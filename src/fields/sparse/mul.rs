use super::num::Zero;
use crate::fields::sparse::basis::{convert_to_base, try_reduce};
use crate::fields::sparse::*;
use crate::fields::{CyclotomicFieldElement, MultiplicativeGroupElement, Q};
use galois::apply_automorphism;
use std::convert::TryInto;
use crate::fields::util::Sign;
use crate::fields::exponent::Exponent;
use crate::fields::rational::Rational;

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
                let new_exp = ((exp1.clone() + exp2.clone()).into(): E) % z1.order.clone();
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
