use crate::fields::sparse_vec::{Number, ExpCoeffMap};
use crate::fields::{AdditiveGroupElement, CyclotomicFieldElement, Q, Z};
use crate::fields::sparse_vec::basis::try_reduce;

impl AdditiveGroupElement for Number {
    fn add(&mut self, rhs: &mut Self) -> &mut Self {
        let mut z1 = self;
        let mut z2 = rhs;
        Self::match_orders(&mut z1, &mut z2);

        let mut i1 = 0;
        let mut i2 = 0;

        let mut result_coeffs = vec![];

        while i1 < z1.coeffs.len() || i2 < z2.coeffs.len() {
            let mut exp1 = &z1.coeffs[i1].0;
            let mut coeff1 = &mut z1.coeffs[i1].1;
            let mut exp2 = &z2.coeffs[i2].0;
            let mut coeff2 = &z2.coeffs[i2].1;

            // the terms of z1 with no matching term in z2, they don't change
            while exp1 < exp2 && i1 < z1.coeffs.len() {
                result_coeffs.push((exp1, coeff1));

                i1 += 1;
                exp1 = &z1.coeffs[i1].0;
                coeff1 = &mut z1.coeffs[i1].1;
            }

            // the terms of z2 with no matching term in z1, these do appear
            // in the result
        }

        self
    }

    fn add_invert(&mut self) -> &mut Self {
        let minus_one = Q::new(Z::from(-1), Z::from(1));
        self.scalar_mul(&minus_one)
    }
}
