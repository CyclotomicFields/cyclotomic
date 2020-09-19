use crate::fields::{FieldElement, AdditiveGroupElement, MultiplicativeGroupElement};
use std::ops::{AddAssign, MulAssign};

pub trait Rational: FieldElement + Clone + From<(i64, u64)> {
    fn zero() -> Self {
        Self::from((0, 1))
    }
    fn is_zero(&self) -> bool;
}

impl<T: Rational> From<i64> for T {
    fn from(n: i64) -> Self {
        Self::from((n, 1))
    }
}

impl AdditiveGroupElement for rug::Rational {
    fn add(&mut self, z: &mut Self) -> &mut Self {
        self.add_assign(z);
        self
    }

    fn add_invert(&mut self) -> &mut Self {
        self.neg_assign();
        self
    }
}

impl MultiplicativeGroupElement for rug::Rational {
    fn mul(&mut self, z: &mut Self) -> &mut Self {
        self.mul_assign(z);
        self
    }

    fn mul_invert(&mut self) -> &mut Self {
        self.recip_mut();
        self
    }
}

impl FieldElement for rug::Rational {
    fn eq(&mut self, other: &mut Self) -> bool {
        self == other
    }
}

impl Rational for rug::Rational {
    fn zero() -> Self {
        rug::Rational::from(0)
    }

    fn is_zero(&self) -> Self {
        self.is_zero()
    }
}

#[derive(Clone, Eq, PartialEq, Debug, Hash)]
pub struct FixedSizeRational {
    pub num: i64,
    denom: u64,
}

/// Not a real rational, but does work if all of your coefficients happen to
/// be representable as 64 bit integers.
impl FixedSizeRational {
    pub fn new(num: i64, denom: u64) -> Self {
        FixedSizeRational {
            num: num,
            denom: denom,
        }
    }

    fn reduce(&mut self) {
        let gcd = num::integer::gcd(self.num, self.denom as i64);
        self.num /= gcd;
        self.denom /= gcd as u64;
    }
}

impl From<(i64, u64)> for FixedSizeRational {
    fn from(p: (i64, u64)) -> Self {
        Rat::new(p.0, p.1)
    }
}

impl AdditiveGroupElement for FixedSizeRational {
    fn add(&mut self, other: &Self) -> &mut Self {
        let common_denom = num::integer::lcm(self.denom, other.denom);
        self.num = self.num * (common_denom / self.denom) as i64
            + other.num * (common_denom / other.denom) as i64;
        self.denom = common_denom;
        self.reduce();
        self
    }

    fn add_invert(&mut self) -> &mut Self {
        self.num *= -1;
        self
    }
}

impl MultiplicativeGroupElement for FixedSizeRational {
    fn mul(&mut self, other: &Self) -> &mut Self {
        self.num *= other.num;
        self.denom *= other.denom;
        self.reduce();
        self
    }

    fn mul_invert(&mut self) -> &mut Self {
        *self = FixedSizeRational::new(
            num::signum(self.num) * self.denom as i64,
            num::abs(self.num) as u64,
        );
        self
    }
}

impl FieldElement for FixedSizeRational {
    fn eq(&mut self, other: &mut Self) -> bool {
        self == other
    }
}

impl Rational for FixedSizeRational {
    fn zero() -> Self {
        FixedSizeRational::from((0, 1))
    }

    fn is_zero(&self) -> Self {
        self.is_zero()
    }
}

