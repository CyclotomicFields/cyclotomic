use crate::fields::sparse_vec::basis::try_reduce;
use crate::fields::sparse_vec::Number;
use crate::fields::{AdditiveGroupElement, CyclotomicFieldElement, Q, Z};

impl AdditiveGroupElement for Number {
    fn add(&mut self, rhs: &mut Self) -> &mut Self {
        let mut z1 = self;
        let mut z2 = rhs;
        Self::match_orders(z1, z2);

        let mut i1 = 0;
        let mut i2 = 0;

        // TODO: what's the best preallocated capacity to use here?
        let mut result_coeffs: Vec<(i64, Q)> = vec![];

        while i1 < z1.coeffs.len() || i2 < z2.coeffs.len() {
            let mut exp1 = z1.coeffs[i1].0;
            let mut coeff1 = &z1.coeffs[i1].1;
            let mut exp2 = z2.coeffs[i2].0;
            let mut coeff2 = &z2.coeffs[i2].1;

            // the terms of z1 with no matching term in z2
            while (exp1 < exp2 || i2 == z2.coeffs.len()) && i1 < z1.coeffs.len() {
                result_coeffs.push((exp1, coeff1.clone()));

                i1 += 1;
                exp1 = z1.coeffs[i1].0;
                coeff1 = &z1.coeffs[i1].1;
            }

            // the terms of z2 with no matching term in z1
            while (exp2 < exp1 || i1 == z1.coeffs.len()) && i2 < z2.coeffs.len() {
                result_coeffs.push((exp2, coeff2.clone()));

                i2 += 1;
                exp2 = z2.coeffs[i2].0;
                coeff2 = &z2.coeffs[i2].1;
            }

            // the matching terms, the ones we have to actually add
            while exp1 == exp2 && i1 < z1.coeffs.len() && i2 < z2.coeffs.len() {
                result_coeffs.push((exp1, coeff1.clone() + coeff2.clone()));

                i1 += 1;
                exp1 = z1.coeffs[i1].0;
                coeff1 = &z1.coeffs[i1].1;

                i2 += 1;
                exp2 = z2.coeffs[i2].0;
                coeff2 = &z2.coeffs[i2].1;
            }
        }

        z1.coeffs = result_coeffs;

        z1
    }

    fn add_invert(&mut self) -> &mut Self {
        let minus_one = Q::new(Z::from(-1), Z::from(1));
        self.scalar_mul(&minus_one)
    }
}
