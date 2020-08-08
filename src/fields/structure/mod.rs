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
    pub basis: Vec<i64>,

    /// \phi(n)
    phi_n: i64,

    /// Let n = \prod_i p_i^{n_i} be a prime factorisation of n. Then factors[i]
    /// = (p_i, n_i).
    factors: Vec<(i64, i64)>,

    zero: Vec<Q>,

    one: Vec<Q>,
}

pub fn write_dense_in_basis(dense: &mut Number, basis: &Vec<i64>) -> Vec<Q> {
    let phi_n = phi(dense.coeffs.len() as i64);
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
    let phi_n = phi(order) as usize;

    // We have to store n^3 space, but maybe if I vectorise this properly it'll
    // be fine?
    let mut structure_constants = vec![vec![vec![Q::from(0); phi_n]; phi_n]; phi_n];

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

            if *p == 2 {
                for bad_exp in (q / 2)..q - 1 + 1 {
                    if math_mod(&i, &q) == math_mod(&((order / q) * bad_exp), &q) {
                        return false;
                    }
                }
            } else {
                for bad_exp in -(q / *p - 1) / 2..(q / *p - 1) / 2 + 1 {
                    if math_mod(&i, &q) == math_mod(&((order / q) * bad_exp), &q) {
                        return false;
                    }
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
        write_dense_in_basis(&mut Number::e(self.order, k), &self.basis)
    }
}

#[derive(Clone, Debug)]
struct SmallOrder(i64);

impl Arbitrary for SmallOrder {
    fn arbitrary<G>(g: &mut G) -> Self
    where
        G: Gen,
    {
        let small_int = g.gen_range(2, 20);
        SmallOrder(small_int)
    }
}

fn random_cyc(field: &CyclotomicField) -> Vec<Q> {
    let mut rng = rand::thread_rng();
    let mut result = vec![];
    for _ in 0..field.phi_n {
        if rng.gen_range(0, 2) == 1 {
            result.push(Q::from(0))
        } else {
            let numerator = rng.gen_range(0, 10);
            let denominator = rng.gen_range(1, 10);
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
        field.basis.len() == phi(small_order.0) as usize
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
}
