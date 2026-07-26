use cyclotomic::fields::sparse::{ExpCoeffMap, Number};
use cyclotomic::fields::{
    AdditiveGroupElement, CyclotomicFieldElement, FieldElement, GenericCyclotomic,
    MultiplicativeGroupElement,
};
use gap_cyclotom_reference::{Context, Cyclotomic};
use rug::Integer;
use std::collections::HashMap;

fn rust_number(order: i64, terms: &[(u32, i64)]) -> Number {
    let mut coefficients = ExpCoeffMap::default();
    for &(exponent, coefficient) in terms {
        coefficients.insert(exponent as i64, rug::Rational::from(coefficient));
    }
    Number::new(&order, &coefficients)
}

fn c_as_rust(value: &Cyclotomic) -> Number {
    let generic = GenericCyclotomic {
        order: Integer::from(value.order()),
        exp_coeffs: value
            .terms()
            .into_iter()
            .map(|(exponent, coefficient)| (Integer::from(exponent), (coefficient, 1_u64)))
            .collect::<HashMap<_, _>>(),
    };
    Number::from_generic(&generic)
}

#[test]
fn products_agree_with_rust_for_seeded_inputs() {
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

            let mut rust_lhs = rust_number(i64::from(order), &lhs_terms);
            let mut rust_rhs = rust_number(i64::from(order), &rhs_terms);
            let rust_product = rust_lhs.mul(&mut rust_rhs);
            assert!(
                c_as_rust(&c_product).eq(rust_product),
                "product mismatch for order {order}, salt {salt}: {c_product:?} vs {rust_product:?}"
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

        let mut rust_left = rust_number(i64::from(left_order), &left_terms);
        let mut rust_right = rust_number(i64::from(right_order), &right_terms);
        let rust_sum = rust_left.add(&mut rust_right);
        assert!(
            c_as_rust(&c_sum).eq(rust_sum),
            "sum mismatch for orders {left_order} and {right_order}"
        );
    }
}
