// Literal Rust translation of the standalone GAP cyclotomic extraction.
//
// GAP is licensed under GPL-2.0-or-later. This derived implementation keeps
// the C extraction's function decomposition, loop order, packed output, dense
// reusable ResultCyc scratch buffer, and checked i64 coefficient arithmetic.
//
// SPDX-License-Identifier: GPL-2.0-or-later

use crate::Error;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Cyclotomic {
    order: u32,
    coefficients: Vec<i64>,
    exponents: Vec<u32>,
}

pub struct Context {
    result_cyc: Vec<i64>,
    last_n: u32,
    phi: u32,
    is_squarefree: bool,
    number_of_primes: u32,
}

impl Default for Context {
    fn default() -> Self {
        Self {
            // GAP initializes ResultCyc with room for 1024 coefficients.
            result_cyc: vec![0; 1024],
            last_n: 0,
            phi: 0,
            is_squarefree: false,
            number_of_primes: 0,
        }
    }
}

fn checked_add(lhs: i64, rhs: i64) -> Result<i64, Error> {
    lhs.checked_add(rhs)
        .ok_or_else(|| Error("coefficient addition overflow".into()))
}

fn checked_sub(lhs: i64, rhs: i64) -> Result<i64, Error> {
    lhs.checked_sub(rhs)
        .ok_or_else(|| Error("coefficient subtraction overflow".into()))
}

fn checked_mul(lhs: i64, rhs: i64) -> Result<i64, Error> {
    lhs.checked_mul(rhs)
        .ok_or_else(|| Error("coefficient multiplication overflow".into()))
}

fn gcd(mut lhs: u32, mut rhs: u32) -> u32 {
    while rhs != 0 {
        let next = lhs % rhs;
        lhs = rhs;
        rhs = next;
    }
    lhs
}

impl Context {
    pub fn new() -> Self {
        Self::default()
    }

    fn grow_result_cyc(&mut self, order: u32) {
        if self.result_cyc.len() < order as usize {
            self.result_cyc.resize(order as usize, 0);
        }
    }

    fn reset_result_cyc(&mut self, order: u32) {
        self.grow_result_cyc(order);
        self.result_cyc[..order as usize].fill(0);
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
                        let coefficient = self.result_cyc[i as usize];
                        if coefficient != 0 {
                            self.result_cyc[i as usize] = 0;
                            for k in 1..p {
                                let exponent = (u64::from(i)
                                    + u64::from(k) * u64::from(n) / u64::from(p))
                                    % u64::from(n);
                                let slot = &mut self.result_cyc[exponent as usize];
                                *slot = checked_sub(*slot, coefficient)?;
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
    fn cyclotomic(&mut self, mut n: u32, hint: u32) -> Result<Cyclotomic, Error> {
        let mut len = 0_usize;
        let mut exponent_gcd = n;
        let mut coefficients_equal = true;
        let mut common_coefficient = None;

        for i in 0..n {
            let coefficient = self.result_cyc[i as usize];
            if coefficient == 0 {
                continue;
            }
            len += 1;
            exponent_gcd = gcd(exponent_gcd, i);
            match common_coefficient {
                None => common_coefficient = Some(coefficient),
                Some(common) if coefficient != common => coefficients_equal = false,
                Some(_) => {}
            }
        }

        if exponent_gcd > 1 {
            let reduced_n = n / exponent_gcd;
            for i in 1..reduced_n {
                self.result_cyc[i as usize] = self.result_cyc[(i * exponent_gcd) as usize];
                self.result_cyc[(i * exponent_gcd) as usize] = 0;
            }
            n = reduced_n;
        }

        let (phi, squarefree, number_of_primes) = self.field_properties(n);
        if len == phi as usize && coefficients_equal && squarefree {
            self.result_cyc[..n as usize].fill(0);
            let mut common = common_coefficient.unwrap_or(0);
            if number_of_primes % 2 != 0 {
                common = checked_sub(0, common)?;
            }
            self.result_cyc[0] = common;
            n = 1;
            len = usize::from(common != 0);
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
                let coefficient = self.result_cyc[((i + n / p) % n) as usize];
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
                let coefficient = self.result_cyc[((i + n / p) % n) as usize];
                self.result_cyc[i as usize] = checked_sub(0, coefficient)?;
                let mut k = i + n / p;
                while k < i + n {
                    self.result_cyc[(k % n) as usize] = 0;
                    k += n / p;
                }
                i += p;
            }
            len /= (p - 1) as usize;

            for i in 1..n / p {
                self.result_cyc[i as usize] = self.result_cyc[(i * p) as usize];
                self.result_cyc[(i * p) as usize] = 0;
            }
            n /= p;
            p += 2;
        }

        let mut coefficients = Vec::with_capacity(len);
        let mut exponents = Vec::with_capacity(len);
        for i in 0..n {
            let coefficient = self.result_cyc[i as usize];
            if coefficient == 0 {
                continue;
            }
            coefficients.push(coefficient);
            exponents.push(i);
            self.result_cyc[i as usize] = 0;
        }
        Ok(Cyclotomic {
            order: n,
            coefficients,
            exponents,
        })
    }

    fn find_common_field(&mut self, nl: u32, nr: u32) -> Result<(u32, u32, u32), Error> {
        let common = u64::from(nl) * u64::from(nr / gcd(nl, nr));
        let n = u32::try_from(common)
            .map_err(|_| Error("common cyclotomic field exceeds uint32_t".into()))?;
        self.grow_result_cyc(n);
        Ok((n / nl, n / nr, n))
    }

    pub fn from_terms(&mut self, order: u32, terms: &[(u32, i64)]) -> Result<Cyclotomic, Error> {
        if order == 0 {
            return Err(Error("cyclotomic order must be positive".into()));
        }
        self.reset_result_cyc(order);
        for &(exponent, coefficient) in terms {
            let slot = &mut self.result_cyc[(exponent % order) as usize];
            *slot = checked_add(*slot, coefficient)?;
        }
        self.convert_to_base(order)?;
        self.cyclotomic(order, 1)
    }

    pub fn root(&mut self, order: u32, exponent: u32) -> Result<Cyclotomic, Error> {
        self.from_terms(order, &[(exponent, 1)])
    }

    pub fn add(&mut self, lhs: &Cyclotomic, rhs: &Cyclotomic) -> Result<Cyclotomic, Error> {
        let (ml, mr, n) = self.find_common_field(lhs.order, rhs.order)?;
        self.result_cyc[..n as usize].fill(0);

        for (&exponent, &coefficient) in lhs.exponents.iter().zip(&lhs.coefficients) {
            self.result_cyc[(exponent * ml) as usize] = coefficient;
        }
        for (&exponent, &coefficient) in rhs.exponents.iter().zip(&rhs.coefficients) {
            let slot = &mut self.result_cyc[(exponent * mr) as usize];
            *slot = checked_add(*slot, coefficient)?;
        }

        if lhs.order % ml != 0 || rhs.order % mr != 0 {
            self.convert_to_base(n)?;
        }
        self.cyclotomic(n, ml * mr)
    }

    pub fn mul(&mut self, lhs: &Cyclotomic, rhs: &Cyclotomic) -> Result<Cyclotomic, Error> {
        // GAP deliberately uses the operand with fewer packed terms as the
        // right operand, then specializes its coefficient before scanning the
        // longer left operand.
        let (left, right) = if lhs.coefficients.len() < rhs.coefficients.len() {
            (rhs, lhs)
        } else {
            (lhs, rhs)
        };
        let (ml, mr, n) = self.find_common_field(left.order, right.order)?;
        self.result_cyc[..n as usize].fill(0);

        for (&right_exponent, &right_coefficient) in right.exponents.iter().zip(&right.coefficients)
        {
            let offset = u64::from(right_exponent) * u64::from(mr) % u64::from(n);
            for (&left_exponent, &left_coefficient) in left.exponents.iter().zip(&left.coefficients)
            {
                let exponent = (offset + u64::from(left_exponent) * u64::from(ml)) % u64::from(n);
                let slot = &mut self.result_cyc[exponent as usize];
                *slot = match right_coefficient {
                    1 => checked_add(*slot, left_coefficient)?,
                    -1 => checked_sub(*slot, left_coefficient)?,
                    _ => checked_add(*slot, checked_mul(left_coefficient, right_coefficient)?)?,
                };
            }
        }

        self.convert_to_base(n)?;
        self.cyclotomic(n, ml * mr)
    }
}

impl Cyclotomic {
    pub fn order(&self) -> u32 {
        self.order
    }

    pub fn terms(&self) -> Vec<(u32, i64)> {
        self.exponents
            .iter()
            .copied()
            .zip(self.coefficients.iter().copied())
            .collect()
    }
}

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
