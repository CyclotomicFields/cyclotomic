use crate::fields::dense::basis::convert_to_base;
use crate::fields::dense::*;
use crate::fields::util::*;
use crate::fields::*;
use num::Zero;

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
            product = convert_to_base(&mut product);
            for k in 0..phi_n {
                structure_constants[i][j][k] = product.coeffs[basis[k] as usize].clone();
            }
        }
    }

    structure_constants
}

fn zumbroich_basis(order: i64) -> Vec<i64> {
    vec![]
}

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
            basis: basis,
            phi_n: phi(order),
            factors: factorise(order),
        }
    }
}