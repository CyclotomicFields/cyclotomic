use gap_cyclotom_reference::literal::Context as I64Context;
use gap_cyclotom_reference::tagged::Context as TaggedContext;
#[cfg(feature = "libgap")]
use gap_cyclotom_reference::tagged::Cyclotomic as TaggedCyclotomic;
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

#[cfg(feature = "libgap")]
fn assert_tagged_matches_gap(
    actual: &TaggedCyclotomic,
    expected: &gap_cyclotom_reference::libgap::Cyclotomic,
    label: &str,
) {
    let expected_order = expected.order().unwrap();
    assert_eq!(actual.order(), expected_order, "{label}: conductor");

    let mut actual_coefficients = vec![(0_i64, 1_i64); expected_order as usize];
    for (exponent, coefficient) in actual.rational_terms() {
        actual_coefficients[exponent as usize] = (
            coefficient
                .numer()
                .to_i64()
                .unwrap_or_else(|| panic!("{label}: numerator does not fit the libgap bridge")),
            coefficient
                .denom()
                .to_i64()
                .unwrap_or_else(|| panic!("{label}: denominator does not fit the libgap bridge")),
        );
    }
    assert_eq!(
        actual_coefficients,
        expected.coefficients().unwrap(),
        "{label}: canonical coefficients"
    );
}

fn rational_terms(order: u32, salt: u32, count: u32) -> Vec<(u32, (i64, i64))> {
    (0..count)
        .map(|index| {
            let exponent = (index * 37 + salt * 11 + index * index) % order;
            let numerator = i64::from((index * 17 + salt * 13) % 41) - 20;
            let denominator = i64::from((index * 7 + salt * 5) % 11 + 1);
            (exponent, (numerator, denominator))
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

#[test]
fn tagged_in_place_addition_reuses_small_values_and_preserves_overflow_fallback() {
    let mut tagged = TaggedContext::new();
    let mut small = tagged.from_integer_terms(5, &[(1, 2)]).unwrap();
    let increment = tagged.from_integer_terms(5, &[(1, 3)]).unwrap();
    tagged.add_assign(&mut small, &increment).unwrap();
    assert!(small.terms()[0].1.is_small());
    assert_eq!(small.terms()[0].1.numerator(), 5);

    let mut overflowing = tagged.from_integer_terms(1, &[(0, i64::MAX)]).unwrap();
    let one = tagged.from_integer_terms(1, &[(0, 1)]).unwrap();
    tagged.add_assign(&mut overflowing, &one).unwrap();
    assert_eq!(
        overflowing.terms()[0].1.numerator(),
        Integer::from(i64::MAX) + 1
    );
    assert!(!overflowing.terms()[0].1.is_small());
}

#[test]
fn tagged_fused_multiply_add_preserves_overflow_fallback() {
    let mut tagged = TaggedContext::new();
    let mut accumulator = tagged.from_integer_terms(1, &[(0, i64::MAX)]).unwrap();
    let one = tagged.from_integer_terms(1, &[(0, 1)]).unwrap();
    tagged.mul_add_assign(&mut accumulator, &one, &one).unwrap();
    assert_eq!(
        accumulator.terms()[0].1.numerator(),
        Integer::from(i64::MAX) + 1
    );
    assert!(!accumulator.terms()[0].1.is_small());
}

#[test]
fn tagged_optimized_operations_match_composition_and_algebraic_laws() {
    let mut tagged = TaggedContext::new();
    for (case_index, (left_order, right_order, accumulator_order)) in
        [(3_u32, 4_u32, 5_u32), (5, 8, 9), (12, 15, 20), (35, 24, 18)]
            .into_iter()
            .enumerate()
    {
        let salt = case_index as u32 + 3;
        let left = tagged
            .from_terms(left_order, &rational_terms(left_order, salt, 7))
            .unwrap();
        let right = tagged
            .from_terms(right_order, &rational_terms(right_order, salt + 11, 8))
            .unwrap();
        let accumulator = tagged
            .from_terms(
                accumulator_order,
                &rational_terms(accumulator_order, salt + 23, 6),
            )
            .unwrap();

        let expected_sum = tagged.add(&left, &right).unwrap();
        let mut in_place_sum = left.clone();
        tagged.add_assign(&mut in_place_sum, &right).unwrap();
        assert_eq!(in_place_sum, expected_sum, "case {case_index}: add_assign");

        let product = tagged.mul(&left, &right).unwrap();
        let expected_fused = tagged.add(&accumulator, &product).unwrap();
        let mut actual_fused = accumulator.clone();
        tagged
            .mul_add_assign(&mut actual_fused, &left, &right)
            .unwrap();
        assert_eq!(
            actual_fused, expected_fused,
            "case {case_index}: fused multiply-add"
        );

        let conjugate_left = tagged.conjugate(&left).unwrap();
        let conjugate_twice = tagged.conjugate(&conjugate_left).unwrap();
        assert_eq!(
            conjugate_twice, left,
            "case {case_index}: conjugation involution"
        );
        let conjugate_product = tagged.conjugate(&product).unwrap();
        let conjugate_right = tagged.conjugate(&right).unwrap();
        let product_of_conjugates = tagged.mul(&conjugate_left, &conjugate_right).unwrap();
        assert_eq!(
            conjugate_product, product_of_conjugates,
            "case {case_index}: conjugation multiplicativity"
        );

        let scaled_once = tagged.scale_fraction(&left, -6, 5).unwrap();
        let scaled_twice = tagged.scale_fraction(&scaled_once, 5, 3).unwrap();
        let scaled_directly = tagged.scale_integer(&left, -2).unwrap();
        assert_eq!(
            scaled_twice, scaled_directly,
            "case {case_index}: scalar composition"
        );

        let right_plus_accumulator = tagged.add(&right, &accumulator).unwrap();
        let distributed_left = tagged.mul(&left, &right_plus_accumulator).unwrap();
        let left_right = tagged.mul(&left, &right).unwrap();
        let left_accumulator = tagged.mul(&left, &accumulator).unwrap();
        let distributed_right = tagged.add(&left_right, &left_accumulator).unwrap();
        assert_eq!(
            distributed_left, distributed_right,
            "case {case_index}: distributivity"
        );
    }
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
fn tagged_operation_matrix_matches_unmodified_gap_canonically() {
    use gap_cyclotom_reference::libgap::Context as LibgapContext;

    let _guard = LIBGAP_TEST_LOCK.lock().unwrap();
    let gap = LibgapContext::new().unwrap();
    let mut tagged = TaggedContext::new();
    let cases = [
        (1_u32, 1_u32, 4_u32),
        (2, 3, 6),
        (3, 4, 8),
        (5, 8, 10),
        (8, 9, 12),
        (12, 15, 15),
        (16, 18, 18),
        (20, 21, 20),
        (24, 35, 24),
        (45, 56, 28),
        (72, 105, 32),
        (77, 120, 36),
        (252, 275, 40),
        (1009, 5, 7),
        (2520, 7, 9),
    ];

    for (case_index, (left_order, right_order, count)) in cases.into_iter().enumerate() {
        let salt = case_index as u32 + 1;
        let left_terms = rational_terms(left_order, salt, count.min(left_order + 3));
        let right_terms = rational_terms(
            right_order,
            salt + 17,
            count.saturating_sub(1).min(right_order + 3),
        );
        let accumulator_order = if case_index % 2 == 0 {
            left_order
        } else {
            right_order
        };
        let accumulator_terms = rational_terms(
            accumulator_order,
            salt + 31,
            (count / 2 + 1).min(accumulator_order + 2),
        );

        let gap_left = gap.from_terms(left_order, &left_terms).unwrap();
        let gap_right = gap.from_terms(right_order, &right_terms).unwrap();
        let gap_accumulator = gap
            .from_terms(accumulator_order, &accumulator_terms)
            .unwrap();
        let tagged_left = tagged.from_terms(left_order, &left_terms).unwrap();
        let tagged_right = tagged.from_terms(right_order, &right_terms).unwrap();
        let tagged_accumulator = tagged
            .from_terms(accumulator_order, &accumulator_terms)
            .unwrap();
        let prefix = format!("case {case_index}, orders {left_order}/{right_order}");

        assert_tagged_matches_gap(&tagged_left, &gap_left, &format!("{prefix}: left input"));
        assert_tagged_matches_gap(&tagged_right, &gap_right, &format!("{prefix}: right input"));

        let gap_sum = gap.add(&gap_left, &gap_right).unwrap();
        let tagged_sum = tagged.add(&tagged_left, &tagged_right).unwrap();
        assert_tagged_matches_gap(&tagged_sum, &gap_sum, &format!("{prefix}: add"));

        let mut tagged_sum_in_place = tagged_left.clone();
        tagged
            .add_assign(&mut tagged_sum_in_place, &tagged_right)
            .unwrap();
        assert_tagged_matches_gap(
            &tagged_sum_in_place,
            &gap_sum,
            &format!("{prefix}: add_assign"),
        );

        let gap_product = gap.mul(&gap_left, &gap_right).unwrap();
        let tagged_product = tagged.mul(&tagged_left, &tagged_right).unwrap();
        assert_tagged_matches_gap(
            &tagged_product,
            &gap_product,
            &format!("{prefix}: multiply"),
        );

        let gap_fused = gap
            .add(&gap_accumulator, &gap.mul(&gap_left, &gap_right).unwrap())
            .unwrap();
        let mut tagged_fused = tagged_accumulator;
        tagged
            .mul_add_assign(&mut tagged_fused, &tagged_left, &tagged_right)
            .unwrap();
        assert_tagged_matches_gap(
            &tagged_fused,
            &gap_fused,
            &format!("{prefix}: fused multiply-add"),
        );

        let conjugate_terms: Vec<_> = left_terms
            .iter()
            .map(|&(exponent, coefficient)| {
                (
                    (left_order - exponent % left_order) % left_order,
                    coefficient,
                )
            })
            .collect();
        let gap_conjugate = gap.from_terms(left_order, &conjugate_terms).unwrap();
        let tagged_conjugate = tagged.conjugate(&tagged_left).unwrap();
        assert_tagged_matches_gap(
            &tagged_conjugate,
            &gap_conjugate,
            &format!("{prefix}: conjugate"),
        );

        let gap_integer_scalar = gap.from_integer_terms(1, &[(0, -7)]).unwrap();
        let gap_integer_scaled = gap.mul(&gap_left, &gap_integer_scalar).unwrap();
        let tagged_integer_scaled = tagged.scale_integer(&tagged_left, -7).unwrap();
        assert_tagged_matches_gap(
            &tagged_integer_scaled,
            &gap_integer_scaled,
            &format!("{prefix}: integer scale"),
        );

        let gap_fraction_scalar = gap.from_terms(1, &[(0, (5, 7))]).unwrap();
        let gap_fraction_scaled = gap.mul(&gap_left, &gap_fraction_scalar).unwrap();
        let tagged_fraction_scaled = tagged.scale_fraction(&tagged_left, 5, 7).unwrap();
        assert_tagged_matches_gap(
            &tagged_fraction_scaled,
            &gap_fraction_scaled,
            &format!("{prefix}: fraction scale"),
        );
    }
}

#[cfg(feature = "libgap")]
#[test]
fn tagged_zero_and_cancellation_cases_match_unmodified_gap() {
    use gap_cyclotom_reference::libgap::Context as LibgapContext;

    let _guard = LIBGAP_TEST_LOCK.lock().unwrap();
    let gap = LibgapContext::new().unwrap();
    let mut tagged = TaggedContext::new();
    for (case_index, (order, terms)) in [
        (1_u32, vec![]),
        (5, vec![(1, (3, 7)), (1, (-3, 7))]),
        (
            12,
            vec![(0, (1, 2)), (4, (1, 2)), (8, (1, 2)), (0, (-1, 1))],
        ),
        (
            30,
            vec![(1, (2, 3)), (7, (-5, 11)), (1, (-2, 3)), (7, (5, 11))],
        ),
    ]
    .into_iter()
    .enumerate()
    {
        let gap_value = gap.from_terms(order, &terms).unwrap();
        let tagged_value = tagged.from_terms(order, &terms).unwrap();
        assert_tagged_matches_gap(
            &tagged_value,
            &gap_value,
            &format!("cancellation case {case_index}"),
        );
        let gap_square = gap.mul(&gap_value, &gap_value).unwrap();
        let tagged_square = tagged.mul(&tagged_value, &tagged_value).unwrap();
        assert_tagged_matches_gap(
            &tagged_square,
            &gap_square,
            &format!("cancellation square {case_index}"),
        );
    }
}

#[cfg(feature = "libgap")]
#[test]
fn tagged_overflow_promotion_round_trips_against_unmodified_gap() {
    use gap_cyclotom_reference::libgap::Context as LibgapContext;

    let _guard = LIBGAP_TEST_LOCK.lock().unwrap();
    let gap = LibgapContext::new().unwrap();
    let mut tagged = TaggedContext::new();

    let gap_max = gap.from_integer_terms(1, &[(0, i64::MAX)]).unwrap();
    let gap_two = gap.from_integer_terms(1, &[(0, 2)]).unwrap();
    let gap_doubled = gap.mul(&gap_max, &gap_two).unwrap();
    let gap_recovered = gap.quo(&gap_doubled, &gap_two).unwrap();

    let tagged_max = tagged.from_integer_terms(1, &[(0, i64::MAX)]).unwrap();
    let tagged_two = tagged.from_integer_terms(1, &[(0, 2)]).unwrap();
    let tagged_doubled = tagged.mul(&tagged_max, &tagged_two).unwrap();
    assert!(!tagged_doubled.terms()[0].1.is_small());
    let tagged_recovered = tagged.scale_fraction(&tagged_doubled, 1, 2).unwrap();
    assert_tagged_matches_gap(
        &tagged_recovered,
        &gap_recovered,
        "overflowing multiplication divided back into bridge range",
    );

    let gap_one = gap.from_integer_terms(1, &[(0, 1)]).unwrap();
    let gap_incremented = gap.add(&gap_max, &gap_one).unwrap();
    let gap_halved = gap.quo(&gap_incremented, &gap_two).unwrap();
    let tagged_one = tagged.from_integer_terms(1, &[(0, 1)]).unwrap();
    let mut tagged_incremented = tagged_max;
    tagged
        .mul_add_assign(&mut tagged_incremented, &tagged_one, &tagged_one)
        .unwrap();
    assert!(!tagged_incremented.terms()[0].1.is_small());
    let tagged_halved = tagged.scale_fraction(&tagged_incremented, 1, 2).unwrap();
    assert_tagged_matches_gap(
        &tagged_halved,
        &gap_halved,
        "overflowing fused multiply-add divided back into bridge range",
    );
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
