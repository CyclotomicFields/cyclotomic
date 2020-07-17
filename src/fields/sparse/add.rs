use crate::fields::sparse::{Number, ExpCoeffMap};
use crate::fields::{AdditiveGroupElement, CyclotomicFieldElement, Q, Z};
use crate::fields::sparse::basis::try_reduce;

impl AdditiveGroupElement for Number {
    /// Simplest possible - term wise addition using hashing.
    ///
    /// Purposely written so it is obviously symmetric in the parameters, thus
    /// commutative by inspection. Of course, there are tests for that.
    fn add(&mut self, rhs: &mut Self) -> &mut Self {
        let mut z1 = self;
        let mut z2 = rhs;
        Self::match_orders(&mut z1, &mut z2);

        // We will never need to reduce here, you can't add low powers of
        // $\zeta_n$ and get higher powers. Higher powers do not exist.
        let mut coeffs = ExpCoeffMap::default();
        for (exp, coeff) in z1.coeffs.clone().into_iter().chain(z2.coeffs.clone()) {
            match coeffs.get(&exp).clone() {
                Some(existing_coeff) => coeffs.insert(exp, coeff + existing_coeff),
                None => coeffs.insert(exp, coeff),
            };
        }

        let mut result = Number::new(z1.order, &coeffs);
        *z1 = result;
        z1
    }

    fn add_invert(&mut self) -> &mut Self {
        let minus_one = Q::new(Z::from(-1), Z::from(1));
        self.scalar_mul(&minus_one)
    }
}
