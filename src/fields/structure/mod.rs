use crate::fields::dense::basis::convert_to_base;
use crate::fields::dense::*;
use crate::fields::util::*;
use crate::fields::*;
use num::Zero;
use quickcheck::{Arbitrary, Gen};
use rand::Rng;

/// This doesn't really fit the same interface as the rest of the fields module,
/// since we don't just have elements on their own, we have structs representing
/// the field and multiply elements using functions on the field.
pub struct CyclotomicField {
    /// I'll refer to the order of the field as n
    order: i64,

    // TODO: This could probably be some kind of incredible vectorized contiguous
    // memory thing, but vec of vec of vecs is good enough to get the tests
    // passing.
    structure_constants: Vec<Vec<Vec<Q>>>,

    /// The basis we are using for the field is: { \zeta_n^{basis[i]} : 0 \leq
    /// i \leq \phi(n) }
    basis: Vec<i64>,

    /// \phi(n)
    phi_n: i64,

    /// Let n = \prod_i p_i^{n_i} be a prime factorisation of n. Then factors[i]
    /// = (p_i, n_i).
    factors: Vec<(i64, i64)>,

    zero: Vec<Q>,

    one: Vec<Q>
}

fn write_dense_in_basis(dense: &mut Number, basis: &Vec<i64>) -> Vec<Q> {
    let phi_n = phi(dense.coeffs.len() as i64);
    let mut result = vec![Q::zero(); phi_n as usize];
    convert_to_base(dense);

    for i in 0..phi_n {
        result[i as usize] = dense.coeffs[basis[i as usize] as usize].clone();
    }

    result
}

/// The structure constants c_ijk are such that (if b_i is the ith basis element)
/// b_i b_j = \sum_k c_ijk b_k.
fn make_structure_constants(order: i64, basis: &Vec<i64>) -> Vec<Vec<Vec<Q>>> {
    let phi_n = phi(order) as usize;

    // We have to store n^3 space, but maybe if I vectorise this properly it'll
    // be fine?
    let mut structure_constants = vec![vec![vec![Q::zero(); phi_n]; phi_n]; phi_n];

    for i in 0..phi_n {
        for j in 0..phi_n {
            let b_i = Number::e(order, basis[i]);
            let b_j = Number::e(order, basis[j]);
            let mut product: Number = b_i.clone().mul(&mut b_j.clone()).clone();
            structure_constants[i][j] = write_dense_in_basis(&mut product, basis);
        }
    }

    structure_constants
}

fn zumbroich_basis(order: i64) -> Vec<i64> {
    let n_div_powers = factorise(order);

    let is_in_basis = |i| {
        for (p, power) in &n_div_powers {
            // the maximal power of p that divides n
            let q: i64 = p.pow(*power as u32);

            // i is in this set (mod q) iff it is not a basis element
            let start_bad = if *p == 2 { q / 2 } else { -(q / p - 1) / 2 };
            let end_bad = if *p == 2 { q - 1 } else { (q / p - 1) / 2 };

            for bad_exp in start_bad..end_bad + 1 {
                if math_mod(&i, &q) == math_mod(&bad_exp, &q) {
                    return false;
                }
            }
        }
        true
    };

    (0..order).filter(|i| is_in_basis(*i)).collect()
}

// TODO: replace with rob's better version
fn factorise(order: i64) -> Vec<(i64, i64)> {
    let n_divisors: Vec<i64> = divisors::get_divisors(order as u64)
        .into_iter()
        .map(|x| x as i64)
        .collect();

    let mut n_div_powers = count_powers(&order, &n_divisors);

    // if it has no divisors smaller than itself, it's prime
    if n_div_powers.is_empty() {
        n_div_powers.push((order, 1));
    }

    n_div_powers
}

impl CyclotomicField {
    pub fn new(order: i64) -> Self {
        let basis = zumbroich_basis(order);
        CyclotomicField {
            order: order,
            structure_constants: make_structure_constants(order, &basis),
            basis: basis.clone(),
            phi_n: phi(order),
            factors: factorise(order),
            zero: write_dense_in_basis(&mut Number::zero_order(order), &basis.clone()),
            one: write_dense_in_basis(&mut Number::one_order(order), &basis.clone()),
        }
    }

    pub fn add(&self, z1: &Vec<Q>, z2: &Vec<Q>) -> Vec<Q> {
        let mut result = vec![Q::zero(); self.phi_n as usize];
        for i in 0..self.phi_n {
            result[i as usize] = &z1[i as usize] + &z2[i as usize];
        }
        result
    }

    pub fn mul(&self, z1: &Vec<Q>, z2: &Vec<Q>) -> Vec<Q> {
        let mut result = vec![Q::zero(); self.phi_n as usize];

        for i in 0..self.phi_n {
            for j in 0..self.phi_n {
                for k in 0..self.phi_n {
                    result[k as usize] += &z1[i as usize]
                        * &z2[j as usize]
                        * &self.structure_constants[i as usize][j as usize][k as usize];
                }
            }
        }

        result
    }
}

#[derive(Clone, Debug)]
struct SmallOrder(i64);

impl Arbitrary for SmallOrder {
    fn arbitrary<G>(g: &mut G) -> Self
    where
        G: Gen,
    {
        let small_int = g.gen_range(2, 30);
        SmallOrder(small_int)
    }
}

fn random_cyc(field: &CyclotomicField) -> Vec<Q> {
    let mut rng = rand::thread_rng();
    let mut result = vec![];
    for _ in 0..field.phi_n {
        let numerator = rng.gen_range(0, 10);
        let denominator = rng.gen_range(1, 10);
        result.push(Q::new(Z::from(numerator), Z::from(denominator)));
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    // TODO: inverses? never heard of em

    #[quickcheck]
    fn zumbroich_basis_has_phi_n_elems(small_order: SmallOrder) -> bool {
        let field = CyclotomicField::new(small_order.0);
        field.basis.len() == phi(small_order.0) as usize
    }

    #[quickcheck]
    fn zero_is_add_identity(small_order: SmallOrder) -> bool {
        let field = CyclotomicField::new(small_order.0);
        let z1 = random_cyc(&field);
        field.add(&z1, &field.zero) == z1
    }

    #[quickcheck]
    fn add_is_associative(small_order: SmallOrder) -> bool {
        let field = CyclotomicField::new(small_order.0);
        let z1 = random_cyc(&field);
        let z2 = random_cyc(&field);
        let z3 = random_cyc(&field);
        field.add(&field.add(&z1, &z2), &z3) == field.add(&z1, &field.add(&z2, &z3))
    }

    #[quickcheck]
    fn add_is_commutative(small_order: SmallOrder) -> bool {
        let field = CyclotomicField::new(small_order.0);
        let z1 = random_cyc(&field);
        let z2 = random_cyc(&field);
        field.add(&z1, &z2) == field.add(&z2, &z1)
    }

    #[quickcheck]
    fn zero_kills_all(small_order: SmallOrder) -> bool {
        let field = CyclotomicField::new(small_order.0);
        let z1 = random_cyc(&field);
        field.mul(&field.zero, &z1) == field.zero
    }

    #[quickcheck]
    fn one_is_mul_identity(small_order: SmallOrder) -> bool {
        let field = CyclotomicField::new(small_order.0);
        let z1 = random_cyc(&field);
        field.mul(&z1, &field.one) == z1
    }

    #[quickcheck]
    fn mul_is_associative(small_order: SmallOrder) -> bool {
        let field = CyclotomicField::new(small_order.0);
        let z1 = random_cyc(&field);
        let z2 = random_cyc(&field);
        let z3 = random_cyc(&field);
        field.mul(&field.mul(&z1, &z2), &z3) == field.mul(&z1, &field.mul(&z2, &z3))
    }

    #[quickcheck]
    fn mul_is_commutative(small_order: SmallOrder) -> bool {
        let field = CyclotomicField::new(small_order.0);
        let z1 = random_cyc(&field);
        let z2 = random_cyc(&field);
        field.mul(&z1, &z2) == field.mul(&z2, &z1)
    }
}
