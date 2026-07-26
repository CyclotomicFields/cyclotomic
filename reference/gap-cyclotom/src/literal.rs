// Literal Rust translation of the standalone GAP cyclotomic extraction.
//
// GAP is licensed under GPL-2.0-or-later. This derived implementation keeps
// the C extraction's function decomposition, loop order, packed output, dense
// reusable ResultCyc scratch buffer, and checked i64 coefficient arithmetic.
//
// SPDX-License-Identifier: GPL-2.0-or-later

use crate::Error;
use rug::Rational;
use std::fmt::Debug;

#[derive(Clone, Debug, Eq, PartialEq)]
#[doc(hidden)]
pub struct KernelCyclotomic<C> {
    order: u32,
    terms: Vec<(u32, C)>,
}

#[doc(hidden)]
pub struct KernelContext<C> {
    result_cyc: Vec<C>,
    last_n: u32,
    phi: u32,
    is_squarefree: bool,
    number_of_primes: u32,
}

#[doc(hidden)]
pub trait Coefficient: Clone + Debug + Eq {
    fn zero() -> Self;
    fn one() -> Self;
    fn is_zero(&self) -> bool;
    fn is_one(&self) -> bool;
    fn is_minus_one(&self) -> bool;
    fn add(&self, rhs: &Self) -> Result<Self, Error>;
    fn sub(&self, rhs: &Self) -> Result<Self, Error>;
    fn mul(&self, rhs: &Self) -> Result<Self, Error>;
    fn add_assign(&mut self, rhs: &Self) -> Result<(), Error> {
        *self = self.add(rhs)?;
        Ok(())
    }
    fn sub_assign(&mut self, rhs: &Self) -> Result<(), Error> {
        *self = self.sub(rhs)?;
        Ok(())
    }
    fn neg(&self) -> Result<Self, Error> {
        Self::zero().sub(self)
    }
}

impl Coefficient for i64 {
    fn zero() -> Self {
        0
    }

    fn one() -> Self {
        1
    }

    fn is_zero(&self) -> bool {
        *self == 0
    }

    fn is_one(&self) -> bool {
        *self == 1
    }

    fn is_minus_one(&self) -> bool {
        *self == -1
    }

    fn add(&self, rhs: &Self) -> Result<Self, Error> {
        self.checked_add(*rhs)
            .ok_or_else(|| Error("coefficient addition overflow".into()))
    }

    fn sub(&self, rhs: &Self) -> Result<Self, Error> {
        self.checked_sub(*rhs)
            .ok_or_else(|| Error("coefficient subtraction overflow".into()))
    }

    fn mul(&self, rhs: &Self) -> Result<Self, Error> {
        self.checked_mul(*rhs)
            .ok_or_else(|| Error("coefficient multiplication overflow".into()))
    }
}

impl Coefficient for Rational {
    fn zero() -> Self {
        Rational::from(0)
    }

    fn one() -> Self {
        Rational::from(1)
    }

    fn is_zero(&self) -> bool {
        self.numer() == &0
    }

    fn is_one(&self) -> bool {
        self == &1
    }

    fn is_minus_one(&self) -> bool {
        self == &-1
    }

    fn add(&self, rhs: &Self) -> Result<Self, Error> {
        let mut value = self.clone();
        value += rhs;
        Ok(value)
    }

    fn sub(&self, rhs: &Self) -> Result<Self, Error> {
        let mut value = self.clone();
        value -= rhs;
        Ok(value)
    }

    fn mul(&self, rhs: &Self) -> Result<Self, Error> {
        if self.is_zero() || rhs.is_zero() {
            return Ok(Self::zero());
        }
        if self.is_one() {
            return Ok(rhs.clone());
        }
        if rhs.is_one() {
            return Ok(self.clone());
        }
        let mut value = self.clone();
        value *= rhs;
        Ok(value)
    }

    fn neg(&self) -> Result<Self, Error> {
        Ok(-self.clone())
    }
}

impl<C: Coefficient> Default for KernelContext<C> {
    fn default() -> Self {
        Self {
            // GAP initializes ResultCyc with room for 1024 coefficients.
            result_cyc: vec![C::zero(); 1024],
            last_n: 0,
            phi: 0,
            is_squarefree: false,
            number_of_primes: 0,
        }
    }
}

fn gcd(mut lhs: u32, mut rhs: u32) -> u32 {
    while rhs != 0 {
        let next = lhs % rhs;
        lhs = rhs;
        rhs = next;
    }
    lhs
}

impl<C: Coefficient> KernelContext<C> {
    pub fn new() -> Self {
        Self::default()
    }

    fn grow_result_cyc(&mut self, order: u32) {
        if self.result_cyc.len() < order as usize {
            self.result_cyc.resize(order as usize, C::zero());
        }
    }

    fn reset_result_cyc(&mut self, order: u32) {
        self.grow_result_cyc(order);
        self.result_cyc[..order as usize].fill(C::zero());
    }

    // Direct translation of GAP's ConvertToBase: eliminate every root outside
    // the Zumbroich basis with 1 + e_p + ... + e_p^(p-1) = 0.
    fn convert_to_base(&mut self, n: u32) -> Result<(), Error> {
        let mut nn = n;
        let mut p = 2;
        while p <= nn {
            if nn % p == 0 {
                let mut q = p;
                while (nn / q) % p == 0 {
                    q *= p;
                }
                nn /= q;

                let start_bad = if p == 2 {
                    i64::from(q / 2)
                } else {
                    -(i64::from(q / p) - 1) / 2
                };
                let end_bad = if p == 2 {
                    i64::from(q) - 1
                } else {
                    (i64::from(q / p) - 1) / 2
                };

                for raw in start_bad..=end_bad {
                    let mut residue = (i64::from(n / q) * raw) % i64::from(q);
                    if residue < 0 {
                        residue += i64::from(q);
                    }

                    let mut i = residue as u32;
                    while i < n {
                        let coefficient = self.result_cyc[i as usize].clone();
                        if !coefficient.is_zero() {
                            self.result_cyc[i as usize] = C::zero();
                            for k in 1..p {
                                let exponent = (u64::from(i)
                                    + u64::from(k) * u64::from(n) / u64::from(p))
                                    % u64::from(n);
                                let slot = &mut self.result_cyc[exponent as usize];
                                slot.sub_assign(&coefficient)?;
                            }
                        }
                        i += q;
                    }
                }
            }
            p += if p == 2 { 1 } else { 2 };
        }
        Ok(())
    }

    fn field_properties(&mut self, n: u32) -> (u32, bool, u32) {
        if n == self.last_n {
            return (self.phi, self.is_squarefree, self.number_of_primes);
        }
        self.last_n = n;
        self.phi = n;
        self.is_squarefree = true;
        self.number_of_primes = 0;
        let mut remaining = n;
        let mut p = 2;

        while p <= remaining {
            if remaining % p == 0 {
                self.phi = self.phi / p * (p - 1);
                self.number_of_primes += 1;
                remaining /= p;
                if remaining % p == 0 {
                    self.is_squarefree = false;
                }
                while remaining % p == 0 {
                    remaining /= p;
                }
            }
            p += 1;
        }
        (self.phi, self.is_squarefree, self.number_of_primes)
    }

    // Direct translation of GAP's Cyclotomic packing and minimal-conductor
    // reduction. `hint` contains prime divisors which may not be removed.
    fn cyclotomic(&mut self, mut n: u32, hint: u32) -> Result<KernelCyclotomic<C>, Error> {
        let mut len = 0_usize;
        let mut exponent_gcd = n;
        let mut coefficients_equal = true;
        let mut common_coefficient = None;

        for i in 0..n {
            let coefficient = self.result_cyc[i as usize].clone();
            if coefficient.is_zero() {
                continue;
            }
            len += 1;
            exponent_gcd = gcd(exponent_gcd, i);
            match common_coefficient {
                None => common_coefficient = Some(coefficient),
                Some(ref common) if coefficient != *common => coefficients_equal = false,
                Some(_) => {}
            }
        }

        if exponent_gcd > 1 {
            let reduced_n = n / exponent_gcd;
            for i in 1..reduced_n {
                self.result_cyc[i as usize] = self.result_cyc[(i * exponent_gcd) as usize].clone();
                self.result_cyc[(i * exponent_gcd) as usize] = C::zero();
            }
            n = reduced_n;
        }

        let (phi, squarefree, number_of_primes) = self.field_properties(n);
        if len == phi as usize && coefficients_equal && squarefree {
            self.result_cyc[..n as usize].fill(C::zero());
            let mut common = common_coefficient.unwrap_or_else(C::zero);
            if number_of_primes % 2 != 0 {
                common = common.neg()?;
            }
            self.result_cyc[0] = common.clone();
            n = 1;
            len = usize::from(!common.is_zero());
        }

        let divisor_gcd = gcd(phi, len as u32);
        let mut nn = n;
        let mut p = 3;
        while p <= nn && p - 1 <= divisor_gcd {
            if nn % p != 0 {
                p += 2;
                continue;
            }
            nn /= p;
            while nn % p == 0 {
                nn /= p;
            }

            if u64::from(n) % (u64::from(p) * u64::from(p)) == 0
                || len % (p - 1) as usize != 0
                || hint % p == 0
            {
                p += 2;
                continue;
            }

            let mut equal_classes = true;
            let mut i = 0;
            while i < n && equal_classes {
                let coefficient = self.result_cyc[((i + n / p) % n) as usize].clone();
                let mut k = i + 2 * n / p;
                while k < i + n {
                    if self.result_cyc[(k % n) as usize] != coefficient {
                        equal_classes = false;
                        break;
                    }
                    k += n / p;
                }
                i += p;
            }

            if !equal_classes {
                p += 2;
                continue;
            }

            let mut i = 0;
            while i < n {
                let coefficient = self.result_cyc[((i + n / p) % n) as usize].clone();
                self.result_cyc[i as usize] = coefficient.neg()?;
                let mut k = i + n / p;
                while k < i + n {
                    self.result_cyc[(k % n) as usize] = C::zero();
                    k += n / p;
                }
                i += p;
            }
            len /= (p - 1) as usize;

            for i in 1..n / p {
                self.result_cyc[i as usize] = self.result_cyc[(i * p) as usize].clone();
                self.result_cyc[(i * p) as usize] = C::zero();
            }
            n /= p;
            p += 2;
        }

        let mut terms = Vec::with_capacity(len);
        for i in 0..n {
            let coefficient = self.result_cyc[i as usize].clone();
            if coefficient.is_zero() {
                continue;
            }
            terms.push((i, coefficient));
            self.result_cyc[i as usize] = C::zero();
        }
        Ok(KernelCyclotomic { order: n, terms })
    }

    fn find_common_field(&mut self, nl: u32, nr: u32) -> Result<(u32, u32, u32), Error> {
        let common = u64::from(nl) * u64::from(nr / gcd(nl, nr));
        let n = u32::try_from(common)
            .map_err(|_| Error("common cyclotomic field exceeds uint32_t".into()))?;
        self.grow_result_cyc(n);
        Ok((n / nl, n / nr, n))
    }

    pub fn from_terms(
        &mut self,
        order: u32,
        terms: &[(u32, C)],
    ) -> Result<KernelCyclotomic<C>, Error> {
        if order == 0 {
            return Err(Error("cyclotomic order must be positive".into()));
        }
        self.reset_result_cyc(order);
        for (exponent, coefficient) in terms {
            let slot = &mut self.result_cyc[(exponent % order) as usize];
            slot.add_assign(coefficient)?;
        }
        self.convert_to_base(order)?;
        self.cyclotomic(order, 1)
    }

    pub fn root(&mut self, order: u32, exponent: u32) -> Result<KernelCyclotomic<C>, Error> {
        self.from_terms(order, &[(exponent, C::one())])
    }

    pub fn add(
        &mut self,
        lhs: &KernelCyclotomic<C>,
        rhs: &KernelCyclotomic<C>,
    ) -> Result<KernelCyclotomic<C>, Error> {
        let (ml, mr, n) = self.find_common_field(lhs.order, rhs.order)?;
        self.result_cyc[..n as usize].fill(C::zero());

        for (exponent, coefficient) in &lhs.terms {
            self.result_cyc[(*exponent * ml) as usize] = coefficient.clone();
        }
        for (exponent, coefficient) in &rhs.terms {
            let slot = &mut self.result_cyc[(*exponent * mr) as usize];
            slot.add_assign(coefficient)?;
        }

        if lhs.order % ml != 0 || rhs.order % mr != 0 {
            self.convert_to_base(n)?;
        }
        self.cyclotomic(n, ml * mr)
    }

    pub fn mul(
        &mut self,
        lhs: &KernelCyclotomic<C>,
        rhs: &KernelCyclotomic<C>,
    ) -> Result<KernelCyclotomic<C>, Error> {
        // GAP deliberately uses the operand with fewer packed terms as the
        // right operand, then specializes its coefficient before scanning the
        // longer left operand.
        let (left, right) = if lhs.terms.len() < rhs.terms.len() {
            (rhs, lhs)
        } else {
            (lhs, rhs)
        };
        let (ml, mr, n) = self.find_common_field(left.order, right.order)?;
        self.result_cyc[..n as usize].fill(C::zero());

        for (right_exponent, right_coefficient) in &right.terms {
            let offset = u64::from(*right_exponent) * u64::from(mr) % u64::from(n);
            for (left_exponent, left_coefficient) in &left.terms {
                let exponent = (offset + u64::from(*left_exponent) * u64::from(ml)) % u64::from(n);
                let slot = &mut self.result_cyc[exponent as usize];
                if right_coefficient.is_one() {
                    slot.add_assign(left_coefficient)?;
                } else if right_coefficient.is_minus_one() {
                    slot.sub_assign(left_coefficient)?;
                } else {
                    slot.add_assign(&left_coefficient.mul(right_coefficient)?)?;
                }
            }
        }

        self.convert_to_base(n)?;
        self.cyclotomic(n, ml * mr)
    }

    pub fn conjugate(&mut self, value: &KernelCyclotomic<C>) -> Result<KernelCyclotomic<C>, Error> {
        let n = value.order;
        self.reset_result_cyc(n);
        for (exponent, coefficient) in &value.terms {
            self.result_cyc[((n - *exponent) % n) as usize] = coefficient.clone();
        }
        self.convert_to_base(n)?;
        self.cyclotomic(n, 1)
    }

    pub fn scale(
        &mut self,
        value: &KernelCyclotomic<C>,
        scalar: &C,
    ) -> Result<KernelCyclotomic<C>, Error> {
        let n = value.order;
        self.reset_result_cyc(n);
        for (exponent, coefficient) in &value.terms {
            self.result_cyc[*exponent as usize] = coefficient.mul(scalar)?;
        }
        // Scalar multiplication preserves a canonical basis representation.
        self.cyclotomic(n, n)
    }
}

impl<C: Coefficient> KernelCyclotomic<C> {
    pub fn order(&self) -> u32 {
        self.order
    }

    pub fn terms(&self) -> Vec<(u32, C)> {
        self.terms.clone()
    }

    pub(crate) fn map_coefficients<D: Coefficient>(
        &self,
        mut map: impl FnMut(&C) -> D,
    ) -> KernelCyclotomic<D> {
        KernelCyclotomic {
            order: self.order,
            terms: self
                .terms
                .iter()
                .map(|(exponent, coefficient)| (*exponent, map(coefficient)))
                .collect(),
        }
    }
}

pub type Context = KernelContext<i64>;
pub type Cyclotomic = KernelCyclotomic<i64>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn literal_port_reduces_roots_and_products() {
        let mut context = Context::new();
        let e6 = context.root(6, 1).unwrap();
        assert_eq!(e6.order(), 3);

        let e2 = context.root(5, 2).unwrap();
        let e4 = context.root(5, 4).unwrap();
        let product = context.mul(&e2, &e4).unwrap();
        assert_eq!(product, context.root(5, 1).unwrap());
    }
}
