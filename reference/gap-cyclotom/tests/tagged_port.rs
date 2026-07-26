use gap_cyclotom_reference::literal::Context as I64Context;
use gap_cyclotom_reference::tagged::Context as TaggedContext;
use rug::Integer;

#[cfg(feature = "libgap")]
static LIBGAP_TEST_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

fn integer_terms(order: u32, count: u32, salt: u32) -> Vec<(u32, i64)> {
    (0..count)
        .map(|index| {
            let exponent = (index * 37 + salt * 11) % order;
            let magnitude = i64::from((index * 13 + salt * 3) % 29 + 1);
            let coefficient = if (index + salt) % 2 == 0 {
                magnitude
            } else {
                -magnitude
            };
            (exponent, coefficient)
        })
        .collect()
}

#[test]
fn tagged_integer_kernel_matches_checked_i64_without_promoting() {
    let cases = [(5_u32, 5_u32), (12, 12), (97, 1), (120, 120), (1009, 11)];
    let mut checked = I64Context::new();
    let mut tagged = TaggedContext::new();

    for (order, count) in cases {
        let left_terms = integer_terms(order, count, 1);
        let right_terms = integer_terms(order, count, 2);
        let checked_left = checked.from_terms(order, &left_terms).unwrap();
        let checked_right = checked.from_terms(order, &right_terms).unwrap();
        let tagged_left = tagged.from_integer_terms(order, &left_terms).unwrap();
        let tagged_right = tagged.from_integer_terms(order, &right_terms).unwrap();

        let checked_product = checked.mul(&checked_left, &checked_right).unwrap();
        let tagged_product = tagged.mul(&tagged_left, &tagged_right).unwrap();
        let tagged_terms: Vec<_> = tagged_product
            .terms()
            .into_iter()
            .map(|(exponent, coefficient)| {
                assert!(coefficient.is_small());
                (exponent, coefficient.numerator().to_i64().unwrap())
            })
            .collect();
        assert_eq!(tagged_product.order(), checked_product.order());
        assert_eq!(tagged_terms, checked_product.terms());
    }
}

#[test]
fn tagged_kernel_promotes_overflowing_cyclotomic_coefficients() {
    let mut tagged = TaggedContext::new();
    let left = tagged.from_integer_terms(1, &[(0, i64::MAX)]).unwrap();
    let right = tagged.from_integer_terms(1, &[(0, 2)]).unwrap();
    let product = tagged.mul(&left, &right).unwrap();
    let terms = product.terms();
    assert_eq!(terms.len(), 1);
    assert!(!terms[0].1.is_small());
    assert_eq!(
        terms[0].1.numerator(),
        Integer::from(i64::MAX) * Integer::from(2)
    );
    assert_eq!(terms[0].1.denominator(), 1);
}

#[test]
fn tagged_fraction_scaling_promotes_without_changing_the_packed_value() {
    let mut tagged = TaggedContext::new();
    let value = tagged.from_integer_terms(5, &[(1, 6), (2, -4)]).unwrap();
    let scaled = tagged.scale_fraction(&value, 1, 2).unwrap();
    let terms: Vec<_> = scaled
        .rational_terms()
        .into_iter()
        .map(|(exponent, coefficient)| {
            (
                exponent,
                (
                    coefficient.numer().to_i64().unwrap(),
                    coefficient.denom().to_i64().unwrap(),
                ),
            )
        })
        .collect();
    assert_eq!(scaled.order(), value.order());
    assert_eq!(terms, vec![(1, (3, 1)), (2, (-2, 1))]);
}

#[cfg(feature = "libgap")]
#[test]
fn tagged_rational_products_equal_unmodified_gap() {
    use gap_cyclotom_reference::libgap::Context as LibgapContext;

    let _guard = LIBGAP_TEST_LOCK.lock().unwrap();
    let gap = LibgapContext::new().unwrap();
    let mut tagged = TaggedContext::new();
    for order in [5_u32, 8, 12, 15, 20, 97] {
        let left_terms = [
            (1 % order, (1, 2)),
            (2 % order, (-2, 3)),
            (4 % order, (5, 7)),
        ];
        let right_terms = [(0, (-3, 5)), (3 % order, (7, 11)), (5 % order, (1, 13))];
        let gap_left = gap.from_terms(order, &left_terms).unwrap();
        let gap_right = gap.from_terms(order, &right_terms).unwrap();
        let gap_product = gap.mul(&gap_left, &gap_right).unwrap();

        let tagged_left = tagged.from_terms(order, &left_terms).unwrap();
        let tagged_right = tagged.from_terms(order, &right_terms).unwrap();
        let tagged_product = tagged.mul(&tagged_left, &tagged_right).unwrap();
        let tagged_terms: Vec<_> = tagged_product
            .rational_terms()
            .into_iter()
            .map(|(exponent, coefficient)| {
                (
                    exponent,
                    (
                        coefficient.numer().to_i64().unwrap(),
                        coefficient.denom().to_i64().unwrap(),
                    ),
                )
            })
            .collect();
        let tagged_as_gap = gap
            .from_terms(tagged_product.order(), &tagged_terms)
            .unwrap();
        assert!(
            tagged_as_gap == gap_product,
            "rational product mismatch at order {order}"
        );
    }
}

#[cfg(feature = "libgap")]
#[test]
fn tagged_character_tensor_decompositions_equal_unmodified_gap() {
    use gap_cyclotom_reference::libgap::{CharacterTableCase, Context as LibgapContext};
    use gap_cyclotom_reference::tagged_character::PreparedCharacterTable;

    let _guard = LIBGAP_TEST_LOCK.lock().unwrap();
    let gap = LibgapContext::new().unwrap();
    let mut tagged = TaggedContext::new();
    for case in CharacterTableCase::ALL {
        let table = gap.character_table(case).unwrap();
        let prepared = PreparedCharacterTable::import(&mut tagged, &table).unwrap();

        for lhs in 0..table.rows.len() {
            for rhs in lhs..table.rows.len() {
                let expected = gap
                    .character_tensor_decomposition(case, lhs, rhs, table.rows.len())
                    .unwrap();
                let actual = prepared
                    .tensor_multiplicities(&mut tagged, lhs, rhs)
                    .unwrap();
                assert_eq!(
                    actual,
                    expected,
                    "{} tensor rows {lhs} and {rhs}",
                    case.name()
                );
            }
        }
    }
}
