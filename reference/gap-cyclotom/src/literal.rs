// Literal Rust translation of the standalone GAP cyclotomic extraction.
//
// GAP is licensed under GPL-2.0-or-later. This derived implementation keeps
// the C extraction's function decomposition, loop order, packed output, dense
// reusable ResultCyc scratch buffer, and checked i64 coefficient arithmetic.
//
// SPDX-License-Identifier: GPL-2.0-or-later

use crate::Error;
use rug::{Integer, Rational};
use smallvec::SmallVec;
use std::collections::HashMap;
use std::fmt::Debug;
use std::ops::{Index, IndexMut};
use std::sync::Arc;

#[derive(Clone, Debug, Eq, PartialEq)]
#[doc(hidden)]
pub struct KernelCyclotomic<C> {
    order: u32,
    terms: SmallVec<[(u32, C); 4]>,
}

#[doc(hidden)]
pub struct KernelContext<C> {
    result_cyc: Scratch<C>,
    last_n: u32,
    phi: u32,
    is_squarefree: bool,
    number_of_primes: u32,
    field_properties_cache: HashMap<u32, (u32, bool, u32)>,
    common_field_cache: HashMap<(u32, u32), (u32, u32, u32)>,
    conversion_plan_cache: HashMap<u32, Arc<[ConversionStep]>>,
}

#[derive(Clone, Copy)]
struct ConversionStep {
    p: u32,
    q: u32,
    start_bad: i64,
    end_bad: i64,
}

struct Scratch<C> {
    values: Vec<C>,
    generations: Vec<u32>,
    generation: u32,
    touched: Vec<usize>,
    active_length: usize,
    track_touched: bool,
}

impl<C: Coefficient> Scratch<C> {
    const TRACKING_THRESHOLD: usize = 4096;

    fn new(capacity: usize) -> Self {
        Self {
            values: vec![C::zero(); capacity],
            generations: vec![0; capacity],
            generation: 0,
            touched: Vec::new(),
            active_length: 0,
            track_touched: false,
        }
    }

    fn grow(&mut self, length: usize) {
        if self.values.len() < length {
            self.values.resize(length, C::zero());
            self.generations.resize(length, 0);
        }
    }

    fn reset(&mut self, length: usize) {
        self.grow(length);
        if self.track_touched {
            for &index in &self.touched {
                self.values[index] = C::zero();
            }
        } else {
            self.values[..self.active_length].fill(C::zero());
        }
        self.generation = self.generation.wrapping_add(1);
        if self.generation == 0 {
            self.generations.fill(0);
            self.generation = 1;
        }
        self.touched.clear();
        self.active_length = length;
        self.track_touched = length > Self::TRACKING_THRESHOLD;
    }

    fn clear_touched(&mut self) {
        if self.track_touched {
            for &index in &self.touched {
                self.values[index] = C::zero();
            }
        } else {
            self.values[..self.active_length].fill(C::zero());
        }
    }

    fn is_scalar(&self) -> bool {
        if self.track_touched {
            self.touched
                .iter()
                .all(|&index| index == 0 || self.values[index].is_zero())
        } else {
            self.values[1..self.active_length]
                .iter()
                .all(Coefficient::is_zero)
        }
    }
}

impl<C: Coefficient> Index<usize> for Scratch<C> {
    type Output = C;

    fn index(&self, index: usize) -> &Self::Output {
        &self.values[index]
    }
}

impl<C: Coefficient> IndexMut<usize> for Scratch<C> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        if self.track_touched && self.generations[index] != self.generation {
            self.generations[index] = self.generation;
            self.values[index] = C::zero();
            self.touched.push(index);
        }
        &mut self.values[index]
    }
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
            result_cyc: Scratch::new(1024),
            last_n: 0,
            phi: 0,
            is_squarefree: false,
            number_of_primes: 0,
            field_properties_cache: HashMap::new(),
            common_field_cache: HashMap::new(),
            conversion_plan_cache: HashMap::new(),
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

fn modular_inverse(value: u64, modulus: u64) -> u64 {
    let (mut old_r, mut r) = (value as i64, modulus as i64);
    let (mut old_s, mut s) = (1_i64, 0_i64);
    while r != 0 {
        let quotient = old_r / r;
        (old_r, r) = (r, old_r - quotient * r);
        (old_s, s) = (s, old_s - quotient * s);
    }
    old_s.rem_euclid(modulus as i64) as u64
}

fn signed_residue(value: i64, modulus: u64) -> u64 {
    (i128::from(value).rem_euclid(i128::from(modulus))) as u64
}

impl<C: Coefficient> KernelContext<C> {
    const CACHE_THRESHOLD: u32 = 256;

    pub fn new() -> Self {
        Self::default()
    }

    fn grow_result_cyc(&mut self, order: u32) {
        self.result_cyc.grow(order as usize);
    }

    fn reset_result_cyc(&mut self, order: u32) {
        self.result_cyc.reset(order as usize);
    }

    fn conversion_plan(&mut self, n: u32) -> Arc<[ConversionStep]> {
        if let Some(plan) = self.conversion_plan_cache.get(&n) {
            return Arc::clone(plan);
        }

        let mut steps = Vec::new();
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
                steps.push(ConversionStep {
                    p,
                    q,
                    start_bad,
                    end_bad,
                });
            }
            p += if p == 2 { 1 } else { 2 };
        }
        let plan: Arc<[ConversionStep]> = steps.into();
        self.conversion_plan_cache.insert(n, Arc::clone(&plan));
        plan
    }

    // Direct translation of GAP's ConvertToBase: eliminate every root outside
    // the Zumbroich basis with 1 + e_p + ... + e_p^(p-1) = 0.
    fn convert_to_base(&mut self, n: u32) -> Result<(), Error> {
        if n <= Self::CACHE_THRESHOLD {
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
            return Ok(());
        }

        let plan = self.conversion_plan(n);
        for step in plan.iter() {
            self.apply_conversion_step(n, *step)?;
        }
        Ok(())
    }

    #[inline(always)]
    fn apply_conversion_step(&mut self, n: u32, step: ConversionStep) -> Result<(), Error> {
        for raw in step.start_bad..=step.end_bad {
            let mut residue = (i64::from(n / step.q) * raw) % i64::from(step.q);
            if residue < 0 {
                residue += i64::from(step.q);
            }
            let mut i = residue as u32;
            while i < n {
                let coefficient = self.result_cyc[i as usize].clone();
                if !coefficient.is_zero() {
                    self.result_cyc[i as usize] = C::zero();
                    for k in 1..step.p {
                        let exponent = (u64::from(i)
                            + u64::from(k) * u64::from(n) / u64::from(step.p))
                            % u64::from(n);
                        let slot = &mut self.result_cyc[exponent as usize];
                        slot.sub_assign(&coefficient)?;
                    }
                }
                i += step.q;
            }
        }
        Ok(())
    }

    fn field_properties(&mut self, n: u32) -> (u32, bool, u32) {
        if n == self.last_n {
            return (self.phi, self.is_squarefree, self.number_of_primes);
        }
        if n > Self::CACHE_THRESHOLD {
            if let Some(&(phi, is_squarefree, number_of_primes)) =
                self.field_properties_cache.get(&n)
            {
                self.last_n = n;
                self.phi = phi;
                self.is_squarefree = is_squarefree;
                self.number_of_primes = number_of_primes;
                return (phi, is_squarefree, number_of_primes);
            }
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
        let properties = (self.phi, self.is_squarefree, self.number_of_primes);
        if n > Self::CACHE_THRESHOLD {
            self.field_properties_cache.insert(n, properties);
        }
        properties
    }

    // Direct translation of GAP's Cyclotomic packing and minimal-conductor
    // reduction. `hint` contains prime divisors which may not be removed.
    fn cyclotomic_into(
        &mut self,
        mut n: u32,
        hint: u32,
        terms: &mut SmallVec<[(u32, C); 4]>,
    ) -> Result<u32, Error> {
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
            self.result_cyc.clear_touched();
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

        terms.clear();
        terms.reserve(len);
        for i in 0..n {
            let coefficient = self.result_cyc[i as usize].clone();
            if coefficient.is_zero() {
                continue;
            }
            terms.push((i, coefficient));
        }
        Ok(n)
    }

    fn cyclotomic(&mut self, n: u32, hint: u32) -> Result<KernelCyclotomic<C>, Error> {
        let mut terms = SmallVec::new();
        let order = self.cyclotomic_into(n, hint, &mut terms)?;
        Ok(KernelCyclotomic { order, terms })
    }

    fn find_common_field(&mut self, nl: u32, nr: u32) -> Result<(u32, u32, u32), Error> {
        let common = u64::from(nl) * u64::from(nr / gcd(nl, nr));
        let n = u32::try_from(common)
            .map_err(|_| Error("common cyclotomic field exceeds uint32_t".into()))?;
        if n > Self::CACHE_THRESHOLD {
            if let Some(&field) = self.common_field_cache.get(&(nl, nr)) {
                self.grow_result_cyc(field.2);
                return Ok(field);
            }
        }
        self.grow_result_cyc(n);
        let field = (n / nl, n / nr, n);
        if n > Self::CACHE_THRESHOLD {
            self.common_field_cache.insert((nl, nr), field);
        }
        Ok(field)
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

    pub(crate) fn from_basis_terms(
        &mut self,
        order: u32,
        terms: &[(u32, C)],
    ) -> Result<KernelCyclotomic<C>, Error> {
        self.reset_result_cyc(order);
        for (exponent, coefficient) in terms {
            self.result_cyc[*exponent as usize] = coefficient.clone();
        }
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
        let (n, hint) = self.add_to_result(lhs, rhs)?;
        self.cyclotomic(n, hint)
    }

    fn add_to_result(
        &mut self,
        lhs: &KernelCyclotomic<C>,
        rhs: &KernelCyclotomic<C>,
    ) -> Result<(u32, u32), Error> {
        let (ml, mr, n) = self.find_common_field(lhs.order, rhs.order)?;
        self.reset_result_cyc(n);

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
        Ok((n, ml * mr))
    }

    pub fn add_assign(
        &mut self,
        lhs: &mut KernelCyclotomic<C>,
        rhs: &KernelCyclotomic<C>,
    ) -> Result<(), Error> {
        let (n, hint) = self.add_to_result(lhs, rhs)?;
        let order = self.cyclotomic_into(n, hint, &mut lhs.terms)?;
        lhs.order = order;
        Ok(())
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
        self.reset_result_cyc(n);

        for (right_exponent, right_coefficient) in &right.terms {
            let offset = *right_exponent * mr;
            for (left_exponent, left_coefficient) in &left.terms {
                let exponent = offset + *left_exponent * ml;
                let exponent = if exponent >= n {
                    exponent - n
                } else {
                    exponent
                };
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

    pub fn mul_add_assign(
        &mut self,
        accumulator: &mut KernelCyclotomic<C>,
        lhs: &KernelCyclotomic<C>,
        rhs: &KernelCyclotomic<C>,
    ) -> Result<(), Error> {
        let common_product =
            u64::from(lhs.order) * u64::from(rhs.order / gcd(lhs.order, rhs.order));
        let common_product = u32::try_from(common_product)
            .map_err(|_| Error("common cyclotomic field exceeds uint32_t".into()))?;
        let n = u64::from(common_product)
            * u64::from(accumulator.order / gcd(common_product, accumulator.order));
        let n = u32::try_from(n)
            .map_err(|_| Error("common cyclotomic field exceeds uint32_t".into()))?;
        self.reset_result_cyc(n);

        let accumulator_scale = n / accumulator.order;
        for (exponent, coefficient) in &accumulator.terms {
            self.result_cyc[(*exponent * accumulator_scale) as usize] = coefficient.clone();
        }

        let (left, right) = if lhs.terms.len() < rhs.terms.len() {
            (rhs, lhs)
        } else {
            (lhs, rhs)
        };
        let left_scale = n / left.order;
        let right_scale = n / right.order;
        for (right_exponent, right_coefficient) in &right.terms {
            let offset = *right_exponent * right_scale;
            for (left_exponent, left_coefficient) in &left.terms {
                let exponent = offset + *left_exponent * left_scale;
                let exponent = if exponent >= n {
                    exponent - n
                } else {
                    exponent
                };
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
        let order = self.cyclotomic_into(n, 1, &mut accumulator.terms)?;
        accumulator.order = order;
        Ok(())
    }

    pub fn sum_products(
        &mut self,
        products: &[(&KernelCyclotomic<C>, &KernelCyclotomic<C>)],
    ) -> Result<KernelCyclotomic<C>, Error> {
        let n = self.sum_products_to_result(products)?;
        self.cyclotomic(n, 1)
    }

    fn sum_products_to_result(
        &mut self,
        products: &[(&KernelCyclotomic<C>, &KernelCyclotomic<C>)],
    ) -> Result<u32, Error> {
        let mut n = 1_u32;
        for (lhs, rhs) in products {
            for order in [lhs.order, rhs.order] {
                let common = u64::from(n) * u64::from(order / gcd(n, order));
                n = u32::try_from(common)
                    .map_err(|_| Error("common cyclotomic field exceeds uint32_t".into()))?;
            }
        }
        self.reset_result_cyc(n);

        for (lhs, rhs) in products {
            let (left, right) = if lhs.terms.len() < rhs.terms.len() {
                (*rhs, *lhs)
            } else {
                (*lhs, *rhs)
            };
            let left_scale = n / left.order;
            let right_scale = n / right.order;
            for (right_exponent, right_coefficient) in &right.terms {
                let offset = *right_exponent * right_scale;
                for (left_exponent, left_coefficient) in &left.terms {
                    let exponent = offset + *left_exponent * left_scale;
                    let exponent = if exponent >= n {
                        exponent - n
                    } else {
                        exponent
                    };
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
        }

        self.convert_to_base(n)?;
        Ok(n)
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

impl KernelContext<i64> {
    const CRT_PRIMES: [u32; 6] = [
        167_772_161,
        469_762_049,
        754_974_721,
        998_244_353,
        1_004_535_809,
        2_013_265_921,
    ];

    pub(crate) fn sum_products_crt_basis(
        &mut self,
        products: &[(&KernelCyclotomic<i64>, &KernelCyclotomic<i64>)],
    ) -> Result<Option<(u32, Vec<Integer>)>, Error> {
        let mut n = 1_u32;
        let mut bound = Integer::from(0);
        for (lhs, rhs) in products {
            for order in [lhs.order, rhs.order] {
                let common = u64::from(n) * u64::from(order / gcd(n, order));
                n = u32::try_from(common)
                    .map_err(|_| Error("common cyclotomic field exceeds uint32_t".into()))?;
            }
            let left_l1 = lhs.terms.iter().fold(Integer::from(0), |sum, (_, value)| {
                sum + Integer::from(value.unsigned_abs())
            });
            let right_l1 = rhs.terms.iter().fold(Integer::from(0), |sum, (_, value)| {
                sum + Integer::from(value.unsigned_abs())
            });
            bound += left_l1 * right_l1;
        }

        let mut remaining = n;
        let mut radical = Integer::from(1);
        let mut factor = 2_u32;
        while factor <= remaining {
            if remaining % factor == 0 {
                radical *= factor;
                while remaining % factor == 0 {
                    remaining /= factor;
                }
            }
            factor += if factor == 2 { 1 } else { 2 };
        }
        bound *= radical;
        bound *= 2;

        if bound == 0 {
            return Ok(Some((n, vec![Integer::from(0); n as usize])));
        }
        let mut combined_modulus = Integer::from(1);
        let mut prime_count = 0;
        for prime in Self::CRT_PRIMES {
            combined_modulus *= prime;
            prime_count += 1;
            if combined_modulus > bound {
                break;
            }
        }
        if combined_modulus <= bound {
            return Ok(None);
        }

        let plan = self.conversion_plan(n);
        let mut residue_sets = Vec::with_capacity(prime_count);
        for &prime in &Self::CRT_PRIMES[..prime_count] {
            let modulus = u64::from(prime);
            let mut coefficients = vec![0_u64; n as usize];
            for (lhs, rhs) in products {
                let left_scale = n / lhs.order;
                let right_scale = n / rhs.order;
                for (left_exponent, left_coefficient) in &lhs.terms {
                    for (right_exponent, right_coefficient) in &rhs.terms {
                        let exponent = *left_exponent * left_scale + *right_exponent * right_scale;
                        let exponent = if exponent >= n {
                            exponent - n
                        } else {
                            exponent
                        };
                        let product = (u128::from(signed_residue(*left_coefficient, modulus))
                            * u128::from(signed_residue(*right_coefficient, modulus)))
                            % u128::from(modulus);
                        coefficients[exponent as usize] =
                            ((u128::from(coefficients[exponent as usize]) + product)
                                % u128::from(modulus)) as u64;
                    }
                }
            }
            for step in plan.iter() {
                for raw in step.start_bad..=step.end_bad {
                    let mut residue = (i64::from(n / step.q) * raw) % i64::from(step.q);
                    if residue < 0 {
                        residue += i64::from(step.q);
                    }
                    let mut i = residue as u32;
                    while i < n {
                        let coefficient = coefficients[i as usize];
                        if coefficient != 0 {
                            coefficients[i as usize] = 0;
                            for k in 1..step.p {
                                let exponent = (u64::from(i)
                                    + u64::from(k) * u64::from(n) / u64::from(step.p))
                                    % u64::from(n);
                                let slot = &mut coefficients[exponent as usize];
                                *slot = (*slot + modulus - coefficient) % modulus;
                            }
                        }
                        i += step.q;
                    }
                }
            }
            residue_sets.push(coefficients);
        }

        let mut prefixes = Vec::with_capacity(prime_count.saturating_sub(1));
        let mut prefix = Integer::from(Self::CRT_PRIMES[0]);
        for &prime in &Self::CRT_PRIMES[1..prime_count] {
            let inverse = modular_inverse(prefix.mod_u(prime) as u64, u64::from(prime));
            prefixes.push((prefix.clone(), prime, inverse));
            prefix *= prime;
        }
        let half_modulus = Integer::from(&combined_modulus >> 1);
        let mut reconstructed = Vec::with_capacity(n as usize);
        for index in 0..n as usize {
            let mut value = Integer::from(residue_sets[0][index]);
            for (stage, (prefix, prime, inverse)) in prefixes.iter().enumerate() {
                let modulus = u64::from(*prime);
                let current = u64::from(value.mod_u(*prime));
                let target = residue_sets[stage + 1][index];
                let delta = ((u128::from((target + modulus - current) % modulus)
                    * u128::from(*inverse))
                    % u128::from(modulus)) as u64;
                value += prefix * delta;
            }
            if value > half_modulus {
                value -= &combined_modulus;
            }
            reconstructed.push(value);
        }
        Ok(Some((n, reconstructed)))
    }

    fn integer_quotient_from_result(&mut self, n: u32, divisor: i64) -> Result<i64, Error> {
        let numerator = self.result_cyc[0];
        if self.result_cyc.is_scalar() {
            if numerator % divisor != 0 {
                return Err(Error("sum of products is not integrally divisible".into()));
            }
            return Ok(numerator / divisor);
        }

        let value = self.cyclotomic(n, 1)?;
        match value.terms.as_slice() {
            [] => Ok(0),
            [(0, numerator)] if value.order == 1 && numerator % divisor == 0 => {
                Ok(numerator / divisor)
            }
            [(0, _)] if value.order == 1 => {
                Err(Error("sum of products is not integrally divisible".into()))
            }
            _ => Err(Error("sum of products is not rational".into())),
        }
    }

    pub fn sum_products_integer_quotient(
        &mut self,
        products: &[(&KernelCyclotomic<i64>, &KernelCyclotomic<i64>)],
        divisor: i64,
    ) -> Result<i64, Error> {
        if divisor == 0 {
            return Err(Error("integer quotient divisor must be nonzero".into()));
        }
        let n = self.sum_products_to_result(products)?;
        self.integer_quotient_from_result(n, divisor)
    }

    pub fn sum_product_rows_integer_quotients(
        &mut self,
        products: &[&KernelCyclotomic<i64>],
        weights: &[&KernelCyclotomic<i64>],
        divisor: i64,
    ) -> Result<Vec<i64>, Error> {
        if divisor == 0 {
            return Err(Error("integer quotient divisor must be nonzero".into()));
        }
        if products.is_empty() || weights.len() % products.len() != 0 {
            return Err(Error("invalid product-row dimensions".into()));
        }

        let mut results = Vec::with_capacity(weights.len() / products.len());
        for row in weights.chunks_exact(products.len()) {
            let mut n = 1_u32;
            for value in products.iter().chain(row.iter()) {
                let common = u64::from(n) * u64::from(value.order / gcd(n, value.order));
                n = u32::try_from(common)
                    .map_err(|_| Error("common cyclotomic field exceeds uint32_t".into()))?;
            }
            self.reset_result_cyc(n);
            for (left, right) in products.iter().zip(row) {
                let (left, right) = if left.terms.len() < right.terms.len() {
                    (*right, *left)
                } else {
                    (*left, *right)
                };
                let left_scale = n / left.order;
                let right_scale = n / right.order;
                for (right_exponent, right_coefficient) in &right.terms {
                    let offset = *right_exponent * right_scale;
                    for (left_exponent, left_coefficient) in &left.terms {
                        let exponent = offset + *left_exponent * left_scale;
                        let exponent = if exponent >= n {
                            exponent - n
                        } else {
                            exponent
                        };
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
            }
            self.convert_to_base(n)?;
            results.push(self.integer_quotient_from_result(n, divisor)?);
        }

        Ok(results)
    }
}

impl<C: Coefficient> KernelCyclotomic<C> {
    pub fn order(&self) -> u32 {
        self.order
    }

    pub fn terms(&self) -> Vec<(u32, C)> {
        self.terms.to_vec()
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
        assert!(!e6.terms.spilled());

        let e2 = context.root(5, 2).unwrap();
        let e4 = context.root(5, 4).unwrap();
        let product = context.mul(&e2, &e4).unwrap();
        assert_eq!(product, context.root(5, 1).unwrap());
    }

    #[test]
    fn fused_multiply_add_matches_separate_operations() {
        let mut context = Context::new();
        for order in [5, 8, 12, 15, 20] {
            let accumulator = context
                .from_terms(order, &[(0, 3), (1, -2), (3, 1)])
                .unwrap();
            let lhs = context.from_terms(order, &[(1, 2), (4, -1)]).unwrap();
            let rhs = context.from_terms(order, &[(0, -3), (2, 1)]).unwrap();
            let product = context.mul(&lhs, &rhs).unwrap();
            let expected = context.add(&accumulator, &product).unwrap();
            let mut actual = accumulator;
            context.mul_add_assign(&mut actual, &lhs, &rhs).unwrap();
            assert_eq!(actual, expected, "order {order}");
        }
    }

    #[test]
    fn delayed_sum_of_products_matches_separate_operations() {
        let mut context = Context::new();
        let values: Vec<_> = [5_u32, 8, 9, 12]
            .into_iter()
            .enumerate()
            .map(|(index, order)| {
                context
                    .from_terms(
                        order,
                        &[
                            (index as u32 % order, index as i64 + 1),
                            ((index as u32 * 3 + 1) % order, index as i64 - 3),
                        ],
                    )
                    .unwrap()
            })
            .collect();
        let pairs = [
            (&values[0], &values[1]),
            (&values[2], &values[3]),
            (&values[0], &values[3]),
        ];
        let mut expected = context.mul(pairs[0].0, pairs[0].1).unwrap();
        for (lhs, rhs) in &pairs[1..] {
            let product = context.mul(lhs, rhs).unwrap();
            expected = context.add(&expected, &product).unwrap();
        }
        let actual = context.sum_products(&pairs).unwrap();
        assert_eq!(actual, expected);
    }

    #[test]
    fn integer_quotient_extracts_scalar_without_packing() {
        let mut context = Context::new();
        let six = context.from_terms(1, &[(0, 6)]).unwrap();
        let two = context.from_terms(1, &[(0, 2)]).unwrap();
        let three = context.from_terms(1, &[(0, 3)]).unwrap();
        let four = context.from_terms(1, &[(0, 4)]).unwrap();
        let pairs = [(&six, &two), (&three, &four)];

        assert_eq!(context.sum_products_integer_quotient(&pairs, 6).unwrap(), 4);
        assert!(context.sum_products_integer_quotient(&pairs, 5).is_err());
    }

    #[test]
    fn generation_stamped_scratch_discards_large_sparse_results_and_errors() {
        let mut context = Context::new();
        context.from_terms(5003, &[(4999, 7)]).unwrap();
        assert!(context
            .from_terms(5003, &[(17, i64::MAX), (17, 1)])
            .is_err());
        let actual = context.from_terms(5003, &[(1, 3)]).unwrap();

        let mut fresh = Context::new();
        let expected = fresh.from_terms(5003, &[(1, 3)]).unwrap();
        assert_eq!(actual, expected);
    }

    #[test]
    fn multi_output_integer_quotients_match_individual_dot_products() {
        let mut context = Context::new();
        let six = context.from_terms(1, &[(0, 6)]).unwrap();
        let three = context.from_terms(1, &[(0, 3)]).unwrap();
        let one = context.from_terms(1, &[(0, 1)]).unwrap();
        let two = context.from_terms(1, &[(0, 2)]).unwrap();
        let four = context.from_terms(1, &[(0, 4)]).unwrap();
        let products = [&six, &three];
        let weights = [&two, &four, &one, &two];

        assert_eq!(
            context
                .sum_product_rows_integer_quotients(&products, &weights, 6)
                .unwrap(),
            vec![4, 2]
        );
    }

    #[test]
    fn crt_basis_reconstruction_matches_exact_overflowing_product() {
        let mut small = Context::new();
        let left = small
            .from_terms(5, &[(0, i64::MAX), (1, i64::MAX - 2)])
            .unwrap();
        let right = small.from_terms(5, &[(0, 2), (2, -3)]).unwrap();
        let (order, coefficients) = small
            .sum_products_crt_basis(&[(&left, &right)])
            .unwrap()
            .unwrap();

        let mut exact = KernelContext::<Rational>::new();
        let basis_terms: Vec<_> = coefficients
            .into_iter()
            .enumerate()
            .filter(|(_, coefficient)| coefficient != &0)
            .map(|(exponent, coefficient)| (exponent as u32, Rational::from(coefficient)))
            .collect();
        let actual = exact.from_basis_terms(order, &basis_terms).unwrap();
        let exact_left = left.map_coefficients(|coefficient| Rational::from(*coefficient));
        let exact_right = right.map_coefficients(|coefficient| Rational::from(*coefficient));
        let expected = exact.sum_products(&[(&exact_left, &exact_right)]).unwrap();
        assert_eq!(actual, expected);
    }
}
