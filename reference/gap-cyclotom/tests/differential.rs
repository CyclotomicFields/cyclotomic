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
