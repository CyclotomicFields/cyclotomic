use crate::fields::dense::basis::convert_to_base;
use crate::fields::dense::*;
use crate::fields::util::*;
use crate::fields::*;
use crate::fields::exponent::Exponent;
use num::Zero;
use quickcheck::{Arbitrary, Gen};
use rand::RngExt;
use std::convert::TryFrom;

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

    /// The same constants as `structure_constants`, but in one contiguous
    /// allocation indexed by `(i * phi_n + j) * phi_n + k`.
    structure_constants_flat: Vec<Q>,

    /// For each basis product `(i, j)`, store only nonzero `(k, c_ijk)` terms.
    structure_constants_sparse: Vec<Vec<(usize, Q)>>,

    /// Integer copy of the flattened structure constants, if all constants fit
    /// in i64.
    structure_constants_i64_flat: Option<Vec<i64>>,

    /// Integer sparse constants grouped by input pair `(i, j)`.
    structure_constants_i64_sparse: Option<Vec<Vec<(usize, i64)>>>,

    /// Integer sparse constants grouped by output coordinate `k`.
    structure_constants_i64_by_output: Option<Vec<Vec<(usize, usize, i64)>>>,

    /// The basis we are using for the field is: { \zeta_n^{basis[i]} : 0 \leq
    /// i \leq \phi(n) }
    pub basis: Vec<i64>,

    /// \phi(n)
    phi_n: i64,

    /// Let n = \prod_i p_i^{n_i} be a prime factorisation of n. Then factors[i]
    /// = (p_i, n_i).
    factors: Vec<(i64, u32)>,

    zero: Vec<Q>,

    one: Vec<Q>,
}

pub struct I64SharedDenElement {
    numerators: Vec<i64>,
    denominator: i128,
}

pub fn write_dense_in_basis(dense: &mut Number, basis: &Vec<i64>) -> Vec<Q> {
    let phi_n = Exponent::phi(&(dense.coeffs.len() as i64));
    let mut result = vec![Q::from(0); phi_n as usize];
    let dense_base = convert_to_base(dense);

    for i in 0..phi_n {
        result[i as usize] = dense_base.coeffs[basis[i as usize] as usize].clone();
    }

    result
}

/// The structure constants c_ijk are such that (if b_i is the ith basis element)
/// b_i b_j = \sum_k c_ijk b_k.
fn make_structure_constants(order: i64, basis: &Vec<i64>) -> Vec<Vec<Vec<Q>>> {
    let phi_n = Exponent::phi(&order) as usize;

    // We have to store n^3 space, but maybe if I vectorise this properly it'll
    // be fine?
    let mut structure_constants = vec![vec![vec![Q::from(0); phi_n]; phi_n]; phi_n];

    for i in 0..phi_n {
        for j in 0..phi_n {
            let b_i = Number::e(&order, &basis[i]);
            let b_j = Number::e(&order, &basis[j]);
            let mut product: Number = b_i.clone().mul(&mut b_j.clone()).clone();
            structure_constants[i][j] = write_dense_in_basis(&mut product, basis);
        }
    }

    structure_constants
}

fn flatten_structure_constants(structure_constants: &Vec<Vec<Vec<Q>>>) -> Vec<Q> {
    let phi_n = structure_constants.len();
    let mut result = Vec::with_capacity(phi_n * phi_n * phi_n);

    for i in 0..phi_n {
        for j in 0..phi_n {
            for k in 0..phi_n {
                result.push(structure_constants[i][j][k].clone());
            }
        }
    }

    result
}

fn sparsify_structure_constants(structure_constants: &Vec<Vec<Vec<Q>>>) -> Vec<Vec<(usize, Q)>> {
    let phi_n = structure_constants.len();
    let mut result = vec![vec![]; phi_n * phi_n];

    for i in 0..phi_n {
        for j in 0..phi_n {
            let row = &mut result[i * phi_n + j];
            for k in 0..phi_n {
                let coeff = &structure_constants[i][j][k];
                if *coeff != 0 {
                    row.push((k, coeff.clone()));
                }
            }
        }
    }

    result
}

fn q_to_i64(q: &Q) -> Option<i64> {
    if *q.denom() != 1 {
        return None;
    }
    q.numer().to_i64()
}

fn integer_flatten_structure_constants(structure_constants: &Vec<Vec<Vec<Q>>>) -> Option<Vec<i64>> {
    let phi_n = structure_constants.len();
    let mut result = Vec::with_capacity(phi_n * phi_n * phi_n);

    for i in 0..phi_n {
        for j in 0..phi_n {
            for k in 0..phi_n {
                result.push(q_to_i64(&structure_constants[i][j][k])?);
            }
        }
    }

    Some(result)
}

fn sparsify_i64_by_input_pair(
    flat_constants: &[i64],
    phi_n: usize,
) -> Vec<Vec<(usize, i64)>> {
    let mut result = vec![vec![]; phi_n * phi_n];

    for i in 0..phi_n {
        for j in 0..phi_n {
            let row = &mut result[i * phi_n + j];
            let row_start = (i * phi_n + j) * phi_n;
            for k in 0..phi_n {
                let coeff = flat_constants[row_start + k];
                if coeff != 0 {
                    row.push((k, coeff));
                }
            }
        }
    }

    result
}

fn sparsify_i64_by_output(flat_constants: &[i64], phi_n: usize) -> Vec<Vec<(usize, usize, i64)>> {
    let mut result = vec![vec![]; phi_n];

    for i in 0..phi_n {
        for j in 0..phi_n {
            let row_start = (i * phi_n + j) * phi_n;
            for k in 0..phi_n {
                let coeff = flat_constants[row_start + k];
                if coeff != 0 {
                    result[k].push((i, j, coeff));
                }
            }
        }
    }

    result
}

fn gcd_i128(mut a: i128, mut b: i128) -> i128 {
    a = a.abs();
    b = b.abs();
    while b != 0 {
        let next = a % b;
        a = b;
        b = next;
    }
    a
}

fn lcm_i128(a: i128, b: i128) -> Option<i128> {
    if a == 0 || b == 0 {
        return Some(0);
    }
    (a / gcd_i128(a, b)).checked_mul(b)
}

fn rational_vector_to_shared_den_i64(z: &[Q]) -> Option<(Vec<i64>, i128)> {
    let mut denominator = 1_i128;
    for coeff in z {
        let coeff_den = coeff.denom().to_i128()?;
        denominator = lcm_i128(denominator, coeff_den)?;
    }

    let mut numerators = Vec::with_capacity(z.len());
    for coeff in z {
        let numerator = coeff.numer().to_i128()?;
        let coeff_den = coeff.denom().to_i128()?;
        let scaled = numerator.checked_mul(denominator / coeff_den)?;
        numerators.push(i64::try_from(scaled).ok()?);
    }

    Some((numerators, denominator))
}

fn rational_from_i128(numerator: i128, denominator: i128) -> Q {
    Q::from((Z::from(numerator), Z::from(denominator)))
}

fn i128_accumulators_to_rationals(accumulators: Vec<i128>, denominator: i128) -> Vec<Q> {
    accumulators
        .into_iter()
        .map(|numerator| rational_from_i128(numerator, denominator))
        .collect()
}

fn zumbroich_basis(order: i64) -> Vec<i64> {
    let n_div_powers = Exponent::factorise(&order);

    let is_in_basis = |i| {
        for (p, power) in &n_div_powers {
            // the maximal power of p that divides n
            let q: i64 = p.pow(*power);

            if *p == 2 {
                for bad_exp in (q / 2)..q - 1 + 1 {
                    if Exponent::math_mod(&i, &q) == Exponent::math_mod(&((order / q) * bad_exp), &q) {
                        return false;
                    }
                }
            } else {
                for bad_exp in -(q / *p - 1) / 2..(q / *p - 1) / 2 + 1 {
                    if Exponent::math_mod(&i, &q) == Exponent::math_mod(&((order / q) * bad_exp), &q) {
                        return false;
                    }
                }
            }
        }
        true
    };

    (0..order).filter(|i| is_in_basis(*i)).collect()
}

fn print_rat_vec(v: &Vec<Q>) -> String {
    let mut result = "[".to_string();

    for q in v {
        result += q.to_string().as_str();
        result += ",";
    }

    result + "]"
}

impl CyclotomicField {
    pub fn new(order: i64) -> Self {
        let basis = zumbroich_basis(order);
        let structure_constants = make_structure_constants(order, &basis);
        let structure_constants_flat = flatten_structure_constants(&structure_constants);
        let structure_constants_sparse = sparsify_structure_constants(&structure_constants);
        let structure_constants_i64_flat =
            integer_flatten_structure_constants(&structure_constants);
        let structure_constants_i64_sparse = structure_constants_i64_flat
            .as_ref()
            .map(|constants| sparsify_i64_by_input_pair(constants, basis.len()));
        let structure_constants_i64_by_output = structure_constants_i64_flat
            .as_ref()
            .map(|constants| sparsify_i64_by_output(constants, basis.len()));
        CyclotomicField {
            order: order,
            structure_constants: structure_constants,
            structure_constants_flat: structure_constants_flat,
            structure_constants_sparse: structure_constants_sparse,
            structure_constants_i64_flat: structure_constants_i64_flat,
            structure_constants_i64_sparse: structure_constants_i64_sparse,
            structure_constants_i64_by_output: structure_constants_i64_by_output,
            basis: basis.clone(),
            phi_n: Exponent::phi(&order),
            factors: Exponent::factorise(&order),
            zero: write_dense_in_basis(&mut Number::zero_order(&order), &basis.clone()),
            one: write_dense_in_basis(&mut Number::one_order(&order), &basis.clone()),
        }
    }

    pub fn add(&self, z1: &Vec<Q>, z2: &Vec<Q>) -> Vec<Q> {
        let mut result = vec![Q::from(0); self.phi_n as usize];
        for i in 0..self.phi_n {
            result[i as usize] = (&z1[i as usize] + &z2[i as usize]).into();
        }
        result
    }

    pub fn mul(&self, z1: &Vec<Q>, z2: &Vec<Q>) -> Vec<Q> {
        let mut result = vec![Q::from(0); self.phi_n as usize];

        for i in 0..self.phi_n {
            for j in 0..self.phi_n {
                for k in 0..self.phi_n {
                    let prod1: Q = (&z1[i as usize] * &z2[j as usize]).into();
                    let prod2: Q = (&prod1
                        * &self.structure_constants[i as usize][j as usize][k as usize])
                        .into();
                    result[k as usize] += prod2;
                }
            }
        }

        result
    }

    pub fn mul_flat(&self, z1: &Vec<Q>, z2: &Vec<Q>) -> Vec<Q> {
        let phi_n = self.phi_n as usize;
        let mut result = vec![Q::from(0); phi_n];

        for i in 0..phi_n {
            for j in 0..phi_n {
                let prod1: Q = (&z1[i] * &z2[j]).into();
                let row_start = (i * phi_n + j) * phi_n;
                for k in 0..phi_n {
                    let prod2: Q = (&prod1 * &self.structure_constants_flat[row_start + k]).into();
                    result[k] += prod2;
                }
            }
        }

        result
    }

    pub fn mul_sparse_constants(&self, z1: &Vec<Q>, z2: &Vec<Q>) -> Vec<Q> {
        let phi_n = self.phi_n as usize;
        let mut result = vec![Q::from(0); phi_n];

        for i in 0..phi_n {
            for j in 0..phi_n {
                let row = &self.structure_constants_sparse[i * phi_n + j];
                if row.is_empty() {
                    continue;
                }

                let prod1: Q = (&z1[i] * &z2[j]).into();
                for (k, constant) in row {
                    let prod2: Q = (&prod1 * constant).into();
                    result[*k] += prod2;
                }
            }
        }

        result
    }

    pub fn mul_sparse_constants_skip_zero_inputs(&self, z1: &Vec<Q>, z2: &Vec<Q>) -> Vec<Q> {
        let phi_n = self.phi_n as usize;
        let mut result = vec![Q::from(0); phi_n];

        for i in 0..phi_n {
            if z1[i] == 0 {
                continue;
            }
            for j in 0..phi_n {
                if z2[j] == 0 {
                    continue;
                }

                let row = &self.structure_constants_sparse[i * phi_n + j];
                if row.is_empty() {
                    continue;
                }

                let prod1: Q = (&z1[i] * &z2[j]).into();
                for (k, constant) in row {
                    let prod2: Q = (&prod1 * constant).into();
                    result[*k] += prod2;
                }
            }
        }

        result
    }

    pub fn mul_i64_flat_shared_den(&self, z1: &Vec<Q>, z2: &Vec<Q>) -> Option<Vec<Q>> {
        let constants = self.structure_constants_i64_flat.as_ref()?;
        let phi_n = self.phi_n as usize;
        let (z1_nums, z1_den) = rational_vector_to_shared_den_i64(z1)?;
        let (z2_nums, z2_den) = rational_vector_to_shared_den_i64(z2)?;
        let denominator = z1_den.checked_mul(z2_den)?;
        let mut result = vec![0_i128; phi_n];

        for i in 0..phi_n {
            for j in 0..phi_n {
                let pair_product = (z1_nums[i] as i128).checked_mul(z2_nums[j] as i128)?;
                let row_start = (i * phi_n + j) * phi_n;
                for k in 0..phi_n {
                    let term = pair_product.checked_mul(constants[row_start + k] as i128)?;
                    result[k] = result[k].checked_add(term)?;
                }
            }
        }

        Some(i128_accumulators_to_rationals(result, denominator))
    }

    pub fn mul_i64_sparse_shared_den(&self, z1: &Vec<Q>, z2: &Vec<Q>) -> Option<Vec<Q>> {
        let constants = self.structure_constants_i64_sparse.as_ref()?;
        let phi_n = self.phi_n as usize;
        let (z1_nums, z1_den) = rational_vector_to_shared_den_i64(z1)?;
        let (z2_nums, z2_den) = rational_vector_to_shared_den_i64(z2)?;
        let denominator = z1_den.checked_mul(z2_den)?;
        let mut result = vec![0_i128; phi_n];

        for i in 0..phi_n {
            if z1_nums[i] == 0 {
                continue;
            }
            for j in 0..phi_n {
                if z2_nums[j] == 0 {
                    continue;
                }
                let pair_product = (z1_nums[i] as i128).checked_mul(z2_nums[j] as i128)?;
                for (k, constant) in &constants[i * phi_n + j] {
                    let term = pair_product.checked_mul(*constant as i128)?;
                    result[*k] = result[*k].checked_add(term)?;
                }
            }
        }

        Some(i128_accumulators_to_rationals(result, denominator))
    }

    pub fn mul_i64_output_grouped_shared_den(&self, z1: &Vec<Q>, z2: &Vec<Q>) -> Option<Vec<Q>> {
        let constants = self.structure_constants_i64_by_output.as_ref()?;
        let phi_n = self.phi_n as usize;
        let (z1_nums, z1_den) = rational_vector_to_shared_den_i64(z1)?;
        let (z2_nums, z2_den) = rational_vector_to_shared_den_i64(z2)?;
        let denominator = z1_den.checked_mul(z2_den)?;
        let mut result = vec![0_i128; phi_n];

        for k in 0..phi_n {
            let mut acc = 0_i128;
            for (i, j, constant) in &constants[k] {
                if z1_nums[*i] == 0 || z2_nums[*j] == 0 {
                    continue;
                }
                let pair_product = (z1_nums[*i] as i128).checked_mul(z2_nums[*j] as i128)?;
                let term = pair_product.checked_mul(*constant as i128)?;
                acc = acc.checked_add(term)?;
            }
            result[k] = acc;
        }

        Some(i128_accumulators_to_rationals(result, denominator))
    }

    pub fn to_i64_shared_den_element(&self, z: &Vec<Q>) -> Option<I64SharedDenElement> {
        let (numerators, denominator) = rational_vector_to_shared_den_i64(z)?;
        Some(I64SharedDenElement {
            numerators,
            denominator,
        })
    }

    pub fn mul_i64_sparse_preconverted(
        &self,
        z1: &I64SharedDenElement,
        z2: &I64SharedDenElement,
    ) -> Option<Vec<Q>> {
        let constants = self.structure_constants_i64_sparse.as_ref()?;
        let phi_n = self.phi_n as usize;
        let denominator = z1.denominator.checked_mul(z2.denominator)?;
        let mut result = vec![0_i128; phi_n];

        for i in 0..phi_n {
            if z1.numerators[i] == 0 {
                continue;
            }
            for j in 0..phi_n {
                if z2.numerators[j] == 0 {
                    continue;
                }
                let pair_product =
                    (z1.numerators[i] as i128).checked_mul(z2.numerators[j] as i128)?;
                for (k, constant) in &constants[i * phi_n + j] {
                    let term = pair_product.checked_mul(*constant as i128)?;
                    result[*k] = result[*k].checked_add(term)?;
                }
            }
        }

        Some(i128_accumulators_to_rationals(result, denominator))
    }

    pub fn mul_i64_output_grouped_preconverted(
        &self,
        z1: &I64SharedDenElement,
        z2: &I64SharedDenElement,
    ) -> Option<Vec<Q>> {
        let constants = self.structure_constants_i64_by_output.as_ref()?;
        let phi_n = self.phi_n as usize;
        let denominator = z1.denominator.checked_mul(z2.denominator)?;
        let mut result = vec![0_i128; phi_n];

        for k in 0..phi_n {
            let mut acc = 0_i128;
            for (i, j, constant) in &constants[k] {
                if z1.numerators[*i] == 0 || z2.numerators[*j] == 0 {
                    continue;
                }
                let pair_product =
                    (z1.numerators[*i] as i128).checked_mul(z2.numerators[*j] as i128)?;
                let term = pair_product.checked_mul(*constant as i128)?;
                acc = acc.checked_add(term)?;
            }
            result[k] = acc;
        }

        Some(i128_accumulators_to_rationals(result, denominator))
    }

    pub fn print(&self, z: &Vec<Q>) -> String {
        let mut str_list: Vec<String> = vec![];
        for i in 0..self.phi_n {
            let exp = self.basis[i as usize];
            let coeff = z[i as usize].clone();
            if coeff != 0 {
                str_list.push(String::from(
                    format!("{} * E({})^{}", coeff, self.order, exp).as_str(),
                ))
            }
        }
        "(".to_string() + &str_list.join(" + ") + ")"
    }

    pub fn e(&self, k: i64) -> Vec<Q> {
        write_dense_in_basis(&mut Number::e(&self.order, &k), &self.basis)
    }
}

#[derive(Clone, Debug)]
struct SmallOrder(i64);

impl Arbitrary for SmallOrder {
    fn arbitrary(g: &mut Gen) -> Self {
        let small_int = arbitrary_i64(g, 2, 20);
        SmallOrder(small_int)
    }
}

fn arbitrary_i64(g: &mut Gen, min: i64, max: i64) -> i64 {
    let width = (max - min) as u64;
    min + (u64::arbitrary(g) % width) as i64
}

fn random_cyc(field: &CyclotomicField) -> Vec<Q> {
    let mut rng = rand::rng();
    let mut result = vec![];
    for _ in 0..field.phi_n {
        if rng.random_range(0..2) == 1 {
            result.push(Q::from(0))
        } else {
            let numerator = rng.random_range(0..10);
            let denominator = rng.random_range(1..10);
            result.push(Q::from((numerator, denominator)));
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    // TODO: inverses? never heard of em

    #[test]
    fn constants_are_correct() {
        let order = 4;
        let field = CyclotomicField::new(order);
        for i in 0..field.phi_n {
            for j in 0..field.phi_n {
                let e_i = field.e(field.basis[i as usize]);
                let e_j = field.e(field.basis[j as usize]);
                let prod = field.mul(&e_i, &e_j);
                println!("---");
                println!("order is {}", order);
                println!(
                    "e_i := {};\ne_j := {};\nprod := {};\n",
                    field.print(&e_i),
                    field.print(&e_j),
                    field.print(&prod)
                );
            }
        }
    }

    #[quickcheck]
    fn zumbroich_basis_has_phi_n_elems(small_order: SmallOrder) -> bool {
        let field = CyclotomicField::new(small_order.0);
        field.basis.len() == Exponent::phi(&small_order.0) as usize
    }

    #[test]
    fn zumbroich_basis_is_correct() {
        assert_eq!(zumbroich_basis(2), vec![0]);
        assert_eq!(zumbroich_basis(3), vec![1, 2]);
        assert_eq!(zumbroich_basis(4), vec![0, 1]);
        assert_eq!(zumbroich_basis(5), vec![1, 2, 3, 4]);
        assert_eq!(zumbroich_basis(6), vec![2, 4]);
        assert_eq!(zumbroich_basis(7), vec![1, 2, 3, 4, 5, 6]);
        assert_eq!(zumbroich_basis(8), vec![0, 1, 2, 3]);
        assert_eq!(zumbroich_basis(9), vec![2, 3, 4, 5, 6, 7]);
        assert_eq!(zumbroich_basis(10), vec![2, 4, 6, 8]);
        assert_eq!(zumbroich_basis(11), vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
        assert_eq!(zumbroich_basis(12), vec![4, 7, 8, 11]);
        assert_eq!(
            zumbroich_basis(13),
            vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
        );
        assert_eq!(zumbroich_basis(14), vec![2, 4, 6, 8, 10, 12]);
        assert_eq!(zumbroich_basis(15), vec![1, 2, 4, 7, 8, 11, 13, 14]);
        assert_eq!(zumbroich_basis(16), vec![0, 1, 2, 3, 4, 5, 6, 7]);
        assert_eq!(
            zumbroich_basis(17),
            vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        );
        assert_eq!(zumbroich_basis(18), vec![4, 6, 8, 10, 12, 14]);
        assert_eq!(
            zumbroich_basis(19),
            vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
        );
        assert_eq!(zumbroich_basis(20), vec![1, 4, 8, 9, 12, 13, 16, 17]);
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
        let left = field.mul(&field.mul(&z1, &z2), &z3);
        let right = field.mul(&z1, &field.mul(&z2, &z3));
        left == right
    }

    #[quickcheck]
    fn mul_is_commutative(small_order: SmallOrder) -> bool {
        let field = CyclotomicField::new(small_order.0);
        let z1 = random_cyc(&field);
        let z2 = random_cyc(&field);
        field.mul(&z1, &z2) == field.mul(&z2, &z1)
    }

    #[quickcheck]
    fn mul_distributes_over_add(small_order: SmallOrder) -> bool {
        let field = CyclotomicField::new(small_order.0);
        let z1 = random_cyc(&field);
        let z2 = random_cyc(&field);
        let z3 = random_cyc(&field);
        let left = field.mul(&z1, &field.add(&z2, &z3));
        let right = field.add(&field.mul(&z1, &z2), &field.mul(&z1, &z3));
        left == right
    }

    #[quickcheck]
    fn flat_constants_match_nested_constants(small_order: SmallOrder) -> bool {
        let field = CyclotomicField::new(small_order.0);
        let z1 = random_cyc(&field);
        let z2 = random_cyc(&field);
        field.mul(&z1, &z2) == field.mul_flat(&z1, &z2)
    }

    #[quickcheck]
    fn sparse_constants_match_nested_constants(small_order: SmallOrder) -> bool {
        let field = CyclotomicField::new(small_order.0);
        let z1 = random_cyc(&field);
        let z2 = random_cyc(&field);
        field.mul(&z1, &z2) == field.mul_sparse_constants(&z1, &z2)
    }

    #[quickcheck]
    fn sparse_constants_skip_zero_inputs_match_nested_constants(small_order: SmallOrder) -> bool {
        let field = CyclotomicField::new(small_order.0);
        let z1 = random_cyc(&field);
        let z2 = random_cyc(&field);
        field.mul(&z1, &z2) == field.mul_sparse_constants_skip_zero_inputs(&z1, &z2)
    }

    #[quickcheck]
    fn i64_flat_shared_den_matches_nested_constants(small_order: SmallOrder) -> bool {
        let field = CyclotomicField::new(small_order.0);
        let z1 = random_cyc(&field);
        let z2 = random_cyc(&field);
        field.mul(&z1, &z2) == field.mul_i64_flat_shared_den(&z1, &z2).unwrap()
    }

    #[quickcheck]
    fn i64_sparse_shared_den_matches_nested_constants(small_order: SmallOrder) -> bool {
        let field = CyclotomicField::new(small_order.0);
        let z1 = random_cyc(&field);
        let z2 = random_cyc(&field);
        field.mul(&z1, &z2) == field.mul_i64_sparse_shared_den(&z1, &z2).unwrap()
    }

    #[quickcheck]
    fn i64_output_grouped_shared_den_matches_nested_constants(small_order: SmallOrder) -> bool {
        let field = CyclotomicField::new(small_order.0);
        let z1 = random_cyc(&field);
        let z2 = random_cyc(&field);
        field.mul(&z1, &z2) == field.mul_i64_output_grouped_shared_den(&z1, &z2).unwrap()
    }

    #[quickcheck]
    fn i64_sparse_preconverted_matches_nested_constants(small_order: SmallOrder) -> bool {
        let field = CyclotomicField::new(small_order.0);
        let z1 = random_cyc(&field);
        let z2 = random_cyc(&field);
        let z1_i64 = field.to_i64_shared_den_element(&z1).unwrap();
        let z2_i64 = field.to_i64_shared_den_element(&z2).unwrap();
        field.mul(&z1, &z2) == field.mul_i64_sparse_preconverted(&z1_i64, &z2_i64).unwrap()
    }

    #[quickcheck]
    fn i64_output_grouped_preconverted_matches_nested_constants(small_order: SmallOrder) -> bool {
        let field = CyclotomicField::new(small_order.0);
        let z1 = random_cyc(&field);
        let z2 = random_cyc(&field);
        let z1_i64 = field.to_i64_shared_den_element(&z1).unwrap();
        let z2_i64 = field.to_i64_shared_den_element(&z2).unwrap();
        field.mul(&z1, &z2)
            == field
                .mul_i64_output_grouped_preconverted(&z1_i64, &z2_i64)
                .unwrap()
    }
}
