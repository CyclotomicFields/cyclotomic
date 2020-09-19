use crate::fields::exponent::Exponent;
use crate::fields::rational::Rational;
use crate::fields::CyclotomicFieldElement;
use crate::fields::Q;

/// The standard character inner product:
/// $\langle \chi_1, \chi_2 \rangle = \sum_g \chi_1(g) \overline{\chi_2(g)}$
pub fn inner_product<T: CyclotomicFieldElement<E, Q>, E: Exponent, Q: Rational>(
    sizes: &Vec<i64>,
    chi1: &Vec<T>,
    chi2: &Vec<T>,
) -> T {
    assert_eq!(chi1.len(), chi2.len());
    assert_eq!(chi1.len(), sizes.len());
    let mut sum = T::zero_order(&E::from(1));

    for i in 0..chi1.len() {
        let mut term = chi1[i].clone();
        term.mul(&mut chi2[i].complex_conjugate())
            .scalar_mul(&Q::from((sizes[i], 1)));
        sum.add(&mut term);
    }

    sum
}
