use cyclotomic::fields::dense;
use cyclotomic::fields::sparse::{self, ExpCoeffMap};
use cyclotomic::fields::structure::{write_dense_in_basis, CyclotomicField};
use cyclotomic::fields::{
    AdditiveGroupElement, CyclotomicFieldElement, FieldElement, GenericCyclotomic,
    MultiplicativeGroupElement,
};
use gap_cyclotom_reference::{Context, Cyclotomic};
use rug::Integer;
use std::collections::HashMap;
#[cfg(feature = "libgap")]
use std::sync::Mutex;

#[cfg(feature = "libgap")]
static LIBGAP_TEST_LOCK: Mutex<()> = Mutex::new(());

#[cfg(feature = "libgap")]
fn libgap_generic(value: &gap_cyclotom_reference::libgap::Cyclotomic) -> GenericCyclotomic {
    GenericCyclotomic {
        order: Integer::from(value.order().unwrap()),
        exp_coeffs: value
            .coefficients()
            .unwrap()
            .into_iter()
            .enumerate()
            .filter(|(_, coefficient)| coefficient.0 != 0)
            .map(|(exponent, (numerator, denominator))| {
                (
                    Integer::from(exponent),
                    (numerator, u64::try_from(denominator).unwrap()),
                )
            })
            .collect(),
    }
}

fn sparse_number(order: i64, terms: &[(u32, i64)]) -> sparse::Number {
    let mut coefficients = ExpCoeffMap::default();
    for &(exponent, coefficient) in terms {
        coefficients.insert(exponent as i64, rug::Rational::from(coefficient));
    }
    sparse::Number::new(&order, &coefficients)
}

fn dense_number(order: i64, terms: &[(u32, i64)]) -> dense::Number {
    let mut coefficients = vec![rug::Rational::from(0); order as usize];
    for &(exponent, coefficient) in terms {
        coefficients[exponent as usize] += coefficient;
    }
    dense::Number::new(&order, &coefficients)
}

fn c_generic(value: &Cyclotomic) -> GenericCyclotomic {
    GenericCyclotomic {
        order: Integer::from(value.order()),
        exp_coeffs: value
            .terms()
            .into_iter()
            .map(|(exponent, coefficient)| (Integer::from(exponent), (coefficient, 1_u64)))
            .collect::<HashMap<_, _>>(),
    }
}

fn structure_as_sparse(
    order: i64,
    field: &CyclotomicField,
    coefficients: &[rug::Rational],
) -> sparse::Number {
    let mut terms = ExpCoeffMap::default();
    for (&exponent, coefficient) in field.basis.iter().zip(coefficients) {
        if coefficient != &0 {
            terms.insert(exponent, coefficient.clone());
        }
    }
    sparse::Number::new(&order, &terms)
}

#[test]
fn same_field_operations_agree_across_all_implementations() {
    let orders = [3_u32, 4, 5, 6, 8, 9, 10, 12, 15, 16, 18, 20, 24, 30];
    let mut context = Context::new().unwrap();

    for order in orders {
        for salt in 0..8_u32 {
            let lhs_terms: Vec<_> = (0..order)
                .filter(|exponent| (exponent + salt) % 3 == 0)
                .map(|exponent| (exponent, (i64::from(exponent + salt) % 5) - 2))
                .filter(|term| term.1 != 0)
                .collect();
            let rhs_terms: Vec<_> = (0..order)
                .filter(|exponent| (2 * exponent + salt) % 5 <= 1)
                .map(|exponent| (exponent, (i64::from(2 * exponent + salt) % 7) - 3))
                .filter(|term| term.1 != 0)
                .collect();

            let c_lhs = context.from_terms(order, &lhs_terms).unwrap();
            let c_rhs = context.from_terms(order, &rhs_terms).unwrap();
            let c_product = context.mul(&c_lhs, &c_rhs).unwrap();
            let c_sum = context.add(&c_lhs, &c_rhs).unwrap();
            let c_product_generic = c_generic(&c_product);
            let c_sum_generic = c_generic(&c_sum);

            let mut sparse_lhs = sparse_number(i64::from(order), &lhs_terms);
            let mut sparse_rhs = sparse_number(i64::from(order), &rhs_terms);
            let sparse_product = sparse_lhs.clone().mul(&mut sparse_rhs.clone()).clone();
            let sparse_sum = sparse_lhs.add(&mut sparse_rhs).clone();
            assert!(
                sparse::Number::from_generic(&c_product_generic).eq(&mut sparse_product.clone()),
                "C/sparse product mismatch for order {order}, salt {salt}"
            );
            assert!(
                sparse::Number::from_generic(&c_sum_generic).eq(&mut sparse_sum.clone()),
                "C/sparse sum mismatch for order {order}, salt {salt}"
            );

            let mut dense_lhs = dense_number(i64::from(order), &lhs_terms);
            let mut dense_rhs = dense_number(i64::from(order), &rhs_terms);
            let dense_product = dense_lhs.clone().mul(&mut dense_rhs.clone()).clone();
            let dense_sum = dense_lhs.add(&mut dense_rhs).clone();
            assert!(
                dense::Number::from_generic(&c_product_generic).eq(&mut dense_product.clone()),
                "C/dense product mismatch for order {order}, salt {salt}"
            );
            assert!(
                dense::Number::from_generic(&c_sum_generic).eq(&mut dense_sum.clone()),
                "C/dense sum mismatch for order {order}, salt {salt}"
            );

            let field = CyclotomicField::new(i64::from(order));
            let structure_lhs = write_dense_in_basis(
                &mut dense_number(i64::from(order), &lhs_terms),
                &field.basis,
            );
            let structure_rhs = write_dense_in_basis(
                &mut dense_number(i64::from(order), &rhs_terms),
                &field.basis,
            );
            let structure_product = field.mul(&structure_lhs, &structure_rhs);
            let structure_sum = field.add(&structure_lhs, &structure_rhs);
            assert!(
                sparse::Number::from_generic(&c_product_generic).eq(&mut structure_as_sparse(
                    i64::from(order),
                    &field,
                    &structure_product,
                )),
                "C/structure product mismatch for order {order}, salt {salt}"
            );
            assert!(
                sparse::Number::from_generic(&c_sum_generic).eq(&mut structure_as_sparse(
                    i64::from(order),
                    &field,
                    &structure_sum,
                )),
                "C/structure sum mismatch for order {order}, salt {salt}"
            );
        }
    }
}

#[test]
fn mixed_field_sums_agree_with_rust() {
    let pairs = [(3_u32, 4_u32), (3, 5), (4, 6), (5, 8), (8, 9), (12, 15)];
    let mut context = Context::new().unwrap();

    for (left_order, right_order) in pairs {
        let left_terms = [(1, 2), (2 % left_order, -1)];
        let right_terms = [(1, -3), (3 % right_order, 2)];
        let c_left = context.from_terms(left_order, &left_terms).unwrap();
        let c_right = context.from_terms(right_order, &right_terms).unwrap();
        let c_sum = context.add(&c_left, &c_right).unwrap();

        let mut rust_left = sparse_number(i64::from(left_order), &left_terms);
        let mut rust_right = sparse_number(i64::from(right_order), &right_terms);
        let rust_sum = rust_left.add(&mut rust_right);
        assert!(
            sparse::Number::from_generic(&c_generic(&c_sum)).eq(rust_sum),
            "sum mismatch for orders {left_order} and {right_order}"
        );
    }
}

#[cfg(feature = "libgap")]
#[test]
fn real_gap_rational_arithmetic_agrees_with_rust() {
    use gap_cyclotom_reference::libgap::Context as LibgapContext;

    let _guard = LIBGAP_TEST_LOCK.lock().unwrap();
    let context = LibgapContext::new().unwrap();
    for order in [7_u32, 8, 12, 15, 20] {
        let lhs_terms = [
            (1 % order, (1, 2)),
            (2 % order, (-2, 3)),
            (4 % order, (5, 7)),
        ];
        let rhs_terms = [(0, (-3, 5)), (3 % order, (7, 11)), (5 % order, (1, 13))];
        let gap_lhs = context.from_terms(order, &lhs_terms).unwrap();
        let gap_rhs = context.from_terms(order, &rhs_terms).unwrap();
        let gap_sum = context.add(&gap_lhs, &gap_rhs).unwrap();
        let gap_product = context.mul(&gap_lhs, &gap_rhs).unwrap();

        let lhs_generic = GenericCyclotomic {
            order: Integer::from(order),
            exp_coeffs: lhs_terms
                .into_iter()
                .map(|(exponent, (numerator, denominator))| {
                    (
                        Integer::from(exponent),
                        (numerator, u64::try_from(denominator).unwrap()),
                    )
                })
                .collect(),
        };
        let rhs_generic = GenericCyclotomic {
            order: Integer::from(order),
            exp_coeffs: rhs_terms
                .into_iter()
                .map(|(exponent, (numerator, denominator))| {
                    (
                        Integer::from(exponent),
                        (numerator, u64::try_from(denominator).unwrap()),
                    )
                })
                .collect(),
        };
        let mut rust_lhs: sparse::Number = sparse::Number::from_generic(&lhs_generic);
        let mut rust_rhs: sparse::Number = sparse::Number::from_generic(&rhs_generic);
        let rust_sum = rust_lhs.clone().add(&mut rust_rhs.clone()).clone();
        let rust_product = rust_lhs.mul(&mut rust_rhs).clone();
        let gap_sum_generic = libgap_generic(&gap_sum);
        let gap_product_generic = libgap_generic(&gap_product);

        assert!(
            sparse::Number::from_generic(&gap_sum_generic).eq(&mut rust_sum.clone()),
            "real GAP/Rust rational sum mismatch for order {order}"
        );
        assert!(
            sparse::Number::from_generic(&gap_product_generic).eq(&mut rust_product.clone()),
            "real GAP/Rust rational product mismatch for order {order}"
        );

        let mut dense_lhs: dense::Number = dense::Number::from_generic(&lhs_generic);
        let mut dense_rhs: dense::Number = dense::Number::from_generic(&rhs_generic);
        let dense_sum = dense_lhs.clone().add(&mut dense_rhs.clone()).clone();
        let dense_product = dense_lhs.mul(&mut dense_rhs).clone();
        assert!(
            dense::Number::from_generic(&gap_sum_generic).eq(&mut dense_sum.clone()),
            "real GAP/dense rational sum mismatch for order {order}"
        );
        assert!(
            dense::Number::from_generic(&gap_product_generic).eq(&mut dense_product.clone()),
            "real GAP/dense rational product mismatch for order {order}"
        );

        let field = CyclotomicField::new(i64::from(order));
        let structure_lhs =
            write_dense_in_basis(&mut dense::Number::from_generic(&lhs_generic), &field.basis);
        let structure_rhs =
            write_dense_in_basis(&mut dense::Number::from_generic(&rhs_generic), &field.basis);
        let structure_sum = field.add(&structure_lhs, &structure_rhs);
        let structure_product = field.mul(&structure_lhs, &structure_rhs);
        assert!(
            sparse::Number::from_generic(&gap_sum_generic).eq(&mut structure_as_sparse(
                i64::from(order),
                &field,
                &structure_sum,
            )),
            "real GAP/structure rational sum mismatch for order {order}"
        );
        assert!(
            sparse::Number::from_generic(&gap_product_generic).eq(&mut structure_as_sparse(
                i64::from(order),
                &field,
                &structure_product,
            )),
            "real GAP/structure rational product mismatch for order {order}"
        );
    }
}

#[cfg(feature = "libgap")]
#[test]
fn stressed_real_gap_products_agree_with_rust() {
    use gap_cyclotom_reference::libgap::Context as LibgapContext;

    let _guard = LIBGAP_TEST_LOCK.lock().unwrap();
    let context = LibgapContext::new().unwrap();
    let cases = [
        (97_u32, 1_u32, true),
        (120, 100, true),
        (1009, 1, false),
        (2520, 1, false),
    ];

    for (order, density_percent, include_structure) in cases {
        let count = (order * density_percent).div_ceil(100).max(1);
        let make_terms = |salt: u32| {
            (0..count)
                .map(|index| {
                    let exponent = (index * 37 + salt * 11) % order;
                    let magnitude = 1_000_000_i64 + i64::from((index + salt * 3) % 97 + 1);
                    let numerator = if (index + salt) % 2 == 0 {
                        magnitude
                    } else {
                        -magnitude
                    };
                    let denominator = [2_i64, 3, 5, 7][((index + salt) % 4) as usize];
                    (exponent, (numerator, denominator))
                })
                .collect::<Vec<_>>()
        };
        let left_terms = make_terms(1);
        let right_terms = make_terms(2);
        let generic = |terms: &[(u32, (i64, i64))]| GenericCyclotomic {
            order: Integer::from(order),
            exp_coeffs: terms
                .iter()
                .map(|&(exponent, (numerator, denominator))| {
                    (
                        Integer::from(exponent),
                        (numerator, u64::try_from(denominator).unwrap()),
                    )
                })
                .collect(),
        };
        let left_generic = generic(&left_terms);
        let right_generic = generic(&right_terms);

        let gap_left = context.from_terms(order, &left_terms).unwrap();
        let gap_right = context.from_terms(order, &right_terms).unwrap();
        let gap_product = context.mul(&gap_left, &gap_right).unwrap();
        let gap_product_generic = libgap_generic(&gap_product);

        let mut sparse_left: sparse::Number = sparse::Number::from_generic(&left_generic);
        let mut sparse_right: sparse::Number = sparse::Number::from_generic(&right_generic);
        let sparse_product = sparse_left.mul(&mut sparse_right).clone();
        assert!(
            sparse::Number::from_generic(&gap_product_generic).eq(&mut sparse_product.clone()),
            "real GAP/sparse stress product mismatch for order {order}, density {density_percent}%"
        );

        let mut dense_left: dense::Number = dense::Number::from_generic(&left_generic);
        let mut dense_right: dense::Number = dense::Number::from_generic(&right_generic);
        let dense_product = dense_left.mul(&mut dense_right).clone();
        assert!(
            dense::Number::from_generic(&gap_product_generic).eq(&mut dense_product.clone()),
            "real GAP/dense stress product mismatch for order {order}, density {density_percent}%"
        );

        if include_structure {
            let field = CyclotomicField::new(i64::from(order));
            let structure_left = write_dense_in_basis(
                &mut dense::Number::from_generic(&left_generic),
                &field.basis,
            );
            let structure_right = write_dense_in_basis(
                &mut dense::Number::from_generic(&right_generic),
                &field.basis,
            );
            let structure_product = field.mul(&structure_left, &structure_right);
            assert!(
                sparse::Number::from_generic(&gap_product_generic).eq(&mut structure_as_sparse(
                    i64::from(order),
                    &field,
                    &structure_product,
                )),
                "real GAP/structure stress product mismatch for order {order}, density {density_percent}%"
            );
        }
    }
}
