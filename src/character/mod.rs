/// The sparse i64 implementation seems like the best all-around for performance.
use crate::fields::sparse::Number;
use crate::fields::CyclotomicFieldElement;

/// We don't aim to have a full algebra system with representations of groups,
/// etc. From the point of view of a character, we only care abotu the sizes
/// of the conjugacy classes, since a character takes a single value on each
/// conjugacy class (characters are class functions).
pub struct ConjugacyClasses(Vec<usize>);

/// A map from conjugacy classes to the value the class function takes.
pub struct Character {
    classes: ConjugacyClasses,
    values: Vec<Number>
}

/// The standard character inner product:
/// $\langle \chi_1, \chi_2 \rangle = \sum_g \chi_1(g) \overline{\chi_2(g)}$
pub fn inner_product(chi1: &Character, chi2: &Character) -> Number {
    // TODO: implement this!
    Number::zero_order(3)
}