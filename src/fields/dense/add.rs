use crate::fields::dense::Number;
use crate::fields::{AdditiveGroupElement, CyclotomicFieldElement, Q, Z};

impl AdditiveGroupElement for Number {
    /// Simplest possible - term wise addition using hashing.
    ///
    /// Purposely written so it is obviously symmetric in the parameters, thus
    /// commutative by inspection. Of course, there are tests for that.
    fn add(&mut self, rhs: &mut Self) -> &mut Self {
        let mut z1 = self;
        let mut z2 = rhs;
        Self::match_orders(&mut z1, &mut z2);

        for exp in 0..z1.order {
            z1.coeffs[exp as usize] += z2.coeffs[exp as usize].clone()
        }

        z1
    }

    fn add_invert(&mut self) -> &mut Self {
        let minus_one = Q::new(Z::from(-1), Z::from(1));
        self.scalar_mul(&minus_one)
    }
}
