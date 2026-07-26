use super::num::Zero;
use crate::fields::sparse::basis::{convert_to_base, try_reduce};
use crate::fields::sparse::*;
use crate::fields::{CyclotomicFieldElement, MultiplicativeGroupElement, Q};
use galois::apply_automorphism;
use std::convert::TryFrom;
use crate::fields::util::Sign;
use crate::fields::exponent::Exponent;
use crate::fields::rational::Rational;

/// Reusable dense accumulator for sparse products. Only touched slots are
/// cleared between calls, so callers can amortize allocation across a workload.
pub struct PackedMulScratch<Q: Rational> {
    coefficients: Vec<Q>,
    touched: Vec<usize>,
}

impl<Q: Rational> PackedMulScratch<Q> {
    pub fn new() -> Self {
        Self {
            coefficients: Vec::new(),
            touched: Vec::new(),
        }
    }

    fn prepare(&mut self, order: usize) {
        for index in self.touched.drain(..) {
            self.coefficients[index] = Q::zero();
        }
        if self.coefficients.len() < order {
            self.coefficients.resize_with(order, Q::zero);
        }
    }
}

impl<Q: Rational> Default for PackedMulScratch<Q> {
    fn default() -> Self {
        Self::new()
    }
}

impl<Q: Rational> Number<i64, Q> {
    pub fn mul_packed(&self, rhs: &Self, scratch: &mut PackedMulScratch<Q>) -> Self {
        if self.order == rhs.order {
            return Self::mul_packed_same_order(self, rhs, scratch);
        }

        let mut left = self.clone();
        let mut right = rhs.clone();
        Self::match_orders(&mut left, &mut right);
        Self::mul_packed_same_order(&left, &right, scratch)
    }

    fn mul_packed_same_order(
        left: &Self,
        right: &Self,
        scratch: &mut PackedMulScratch<Q>,
    ) -> Self {
        let order = usize::try_from(left.order).expect("cyclotomic order fits usize");
        scratch.prepare(order);
        for (left_exp, left_coeff) in &left.coeffs {
            for (right_exp, right_coeff) in &right.coeffs {
                let exponent = (left_exp + right_exp).rem_euclid(left.order) as usize;
                let mut product = left_coeff.clone();
                product.mul(&mut right_coeff.clone());
                if scratch.coefficients[exponent].is_zero() {
                    scratch.touched.push(exponent);
                }
                scratch.coefficients[exponent].add(&mut product);
            }
        }

        let mut coefficients = ExpCoeffMap::default();
        for &index in &scratch.touched {
            let coefficient = &scratch.coefficients[index];
            if !coefficient.is_zero() {
                coefficients.insert(index as i64, coefficient.clone());
            }
        }
        Number::new(&left.order, &coefficients)
    }
}

impl<E, Q> MultiplicativeGroupElement for Number<E, Q> where E: Exponent, Q: Rational {
    /// Multiplies term by term, not bothering to do anything interesting.
    fn mul(&mut self, rhs: &mut Self) -> &mut Self {
        let z1 = self;
        let z2 = rhs;
        Self::match_orders(z1, z2);

        let mut result = Number::zero_order(&z1.order);

        // This order is almost certainly not optimal. But you know, whatever.
        // TODO: make it gooder
        result.order = z1.order.clone();
        for (exp1, coeff1) in &mut z1.coeffs {
            for (exp2, coeff2) in &mut z2.coeffs {
                let new_exp = (exp1.clone() + exp2.clone()) % z1.order.clone();
                let new_coeff = coeff1.clone().mul(coeff2).clone();
                add_single(&mut result.coeffs, &new_exp, &new_coeff, Sign::Plus);
            }
        }

        *z1 = result;
        z1
    }

    /// Gives the inverse of $z$ using the product of Galois conjugates.
    ///
    /// I don't think there's a "trivial" or "stupid" way of doing this.
    /// The product of the Galois conjugates is rational, we can normalise
    /// to get the multiplicative inverse.
    fn mul_invert(&mut self) -> &mut Self {
        let z = self.clone();
        // Let $L = \mathbb{Q}(\zeta_n), K = \mathbb{Q}$.
        // Then $L/K$ is a degree $\phi(n)$ extension.
        let n = &z.order;

        // The Galois group $G = \text{Aut}(L/K)$ has order $\phi(n)$. The
        // elements are the automorphisms $\zeta_n \mapsto \zeta_n^i$ for all
        // $1 \leq i \leq n-1$ coprime to $n$.

        // This is the product except for the term for $t = \id_L$.

        let mut x = Number::one_order(n);

        let mut i = E::from(2);
        while &i != n {
            if Exponent::gcd(n ,&i) == E::from(1) {
                x.mul(&mut apply_automorphism(&z, &i));
            }
            i =  i + E::from(1);
        }

        // The full product:
        let mut q_cyc = convert_to_base(z.clone().mul(&mut x.clone()));
        try_reduce(&mut q_cyc);
        println!("q_cyc = {:?}", q_cyc);

        assert_eq!(q_cyc.order, E::from(1));
        let q_rat = q_cyc.coeffs.get_mut(&E::from(0)).unwrap();
        println!("q_rat = {:?}", q_rat);

        if q_rat.is_zero() {
            panic!("can't invert zero!");
        }
        q_rat.mul_invert();

        let z_inv = x.scalar_mul(q_rat);
        *self = z_inv.clone();
        self
    }
}

#[cfg(test)]
mod packed_tests {
    use super::*;
    use crate::fields::rational::HybridRational;

    #[test]
    fn packed_products_match_the_general_implementation_and_reuse_scratch() {
        let left = Number::<i64, HybridRational>::new(
            &5,
            &[(1, 2.into()), (4, (-3).into())]
                .iter()
                .cloned()
                .collect(),
        );
        let right = Number::<i64, HybridRational>::new(
            &3,
            &[(0, 7.into()), (2, 1.into())]
                .iter()
                .cloned()
                .collect(),
        );
        let mut expected = left.clone();
        expected.mul(&mut right.clone());

        let mut scratch = PackedMulScratch::new();
        let mut actual = left.mul_packed(&right, &mut scratch);
        assert!(actual.eq(&mut expected));

        let mut second = right.mul_packed(&right, &mut scratch);
        let mut expected_second = right.clone();
        expected_second.mul(&mut right.clone());
        assert!(second.eq(&mut expected_second));
    }
}
