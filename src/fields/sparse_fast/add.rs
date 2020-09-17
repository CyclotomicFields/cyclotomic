use crate::fields::sparse_fast::{ExpCoeffMap, Number};
use crate::fields::{AdditiveGroupElement, CyclotomicFieldElement, Q, Z};
use crate::fields::exponent::Exponent;
use std::ops::AddAssign;

impl<E> AdditiveGroupElement for Number<E> where E: Exponent {
    /// Simplest possible - term wise addition using hashing.
    ///
    /// Purposely written so it is obviously symmetric in the parameters, thus
    /// commutative by inspection. Of course, there are tests for that.
    fn add(&mut self, rhs: &mut Self) -> &mut Self {
        let mut z1 = self;
        let mut z2 = rhs;
        Self::match_orders(&mut z1, &mut z2);

        for (exp, coeff) in &z2.coeffs {
            match z1.coeffs.get_mut(&exp) {
                Some(existing_coeff) => existing_coeff.add_assign(coeff),
                None => { z1.coeffs.insert(exp.clone(), coeff.clone()); }
            };
        }

        z1
    }

    fn add_invert(&mut self) -> &mut Self {
        let minus_one = Q::from(-1);
        self.scalar_mul(&minus_one)
    }
}
