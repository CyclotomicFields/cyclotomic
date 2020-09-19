use crate::fields::{AdditiveGroupElement, FieldElement, MultiplicativeGroupElement, Z};
use std::fmt;
use std::fmt::{Debug, Display};
use std::hash::{Hash, Hasher};
use std::mem::transmute;
use std::ops::{AddAssign, MulAssign};

pub trait Rational:
    FieldElement + Clone + From<(i64, u64)> + From<i64> + Display + Send + Debug + PartialEq + Eq + Hash
{
    fn zero() -> Self {
        Self::from((0, 1))
    }

    fn is_zero(&self) -> bool;

    fn numer(&self) -> Z;

    fn denom(&self) -> Z;
}

// TODO: there are a lot of clones due to the &mut interface, need to change

impl AdditiveGroupElement for rug::Rational {
    fn add(&mut self, z: &mut Self) -> &mut Self {
        *self += z.clone();
        self
    }

    fn add_invert(&mut self) -> &mut Self {
        *self = self.as_neg().clone();
        self
    }
}

impl MultiplicativeGroupElement for rug::Rational {
    fn mul(&mut self, z: &mut Self) -> &mut Self {
        *self *= z.clone();
        self
    }

    fn mul_invert(&mut self) -> &mut Self {
        self.recip_mut();
        self
    }
}

impl FieldElement for rug::Rational {
    fn eq(&mut self, other: &mut Self) -> bool {
        *self == *other
    }
}

impl Rational for rug::Rational {
    fn zero() -> Self {
        rug::Rational::from(0)
    }

    fn is_zero(&self) -> bool {
        self.is_zero()
    }

    fn numer(&self) -> Z {
        self.numer().clone()
    }

    fn denom(&self) -> Z {
        self.denom().clone()
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

impl From<i64> for FixedSizeRational {
    fn from(n: i64) -> Self {
        Self::from((n, 1))
    }
}

impl Display for FixedSizeRational {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}/{}", self.num, self.denom)
    }
}

impl From<(i64, u64)> for FixedSizeRational {
    fn from(p: (i64, u64)) -> Self {
        FixedSizeRational::new(p.0, p.1)
    }
}

impl AdditiveGroupElement for FixedSizeRational {
    fn add(&mut self, other: &mut Self) -> &mut Self {
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
    fn mul(&mut self, other: &mut Self) -> &mut Self {
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
        *self == *other
    }
}

impl Rational for FixedSizeRational {
    fn zero() -> Self {
        FixedSizeRational::from((0, 1))
    }

    fn is_zero(&self) -> bool {
        self.is_zero()
    }

    fn numer(&self) -> Z {
        Z::from(self.num)
    }

    fn denom(&self) -> Z {
        Z::from(self.denom)
    }
}

// floats aren't rationals, but they are very fast...
#[derive(Debug, Clone, PartialEq)]
pub struct FloatRational(f64);

impl Eq for FloatRational {}

// maybe a hack
impl Hash for FloatRational {
    fn hash<H: Hasher>(&self, state: &mut H) {
        let x: u64 = unsafe { transmute(self.0) };
        x.hash(state);
    }
}

impl Display for FloatRational {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl From<(i64, u64)> for FloatRational {
    fn from(p: (i64, u64)) -> Self {
        FloatRational(p.0 as f64 / p.1 as f64)
    }
}

impl From<i64> for FloatRational {
    fn from(n: i64) -> Self {
        FloatRational(n as f64)
    }
}

impl AdditiveGroupElement for FloatRational {
    fn add(&mut self, other: &mut Self) -> &mut Self {
        self.0 += other.0;
        self
    }

    fn add_invert(&mut self) -> &mut Self {
        self.0 *= -1_f64;
        self
    }
}

impl MultiplicativeGroupElement for FloatRational {
    fn mul(&mut self, other: &mut Self) -> &mut Self {
        self.0 *= other.0;
        self
    }

    fn mul_invert(&mut self) -> &mut Self {
        self.0 = 1_f64 / self.0;
        self
    }
}

impl FieldElement for FloatRational {
    fn eq(&mut self, other: &mut Self) -> bool {
        *self == *other
    }
}

impl Rational for FloatRational {
    fn zero() -> Self {
        FloatRational(0_f64)
    }

    fn is_zero(&self) -> bool {
        self.0 == 0_f64
    }

    // TODO: fix the API so I don't have to lie here
    fn numer(&self) -> Z {
        // this is a hack, possible only because I know this is usually used to
        // get an integer from the rational.
        Z::from(self.0 as i64)
    }

    fn denom(&self) -> Z {
        // this is almost always wrong
        Z::from(1)
    }
}
