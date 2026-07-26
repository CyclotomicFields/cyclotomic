use gap_cyclotom_reference::literal::Context as RustContext;
use gap_cyclotom_reference::Context as CContext;

fn terms(order: u32, count: u32, salt: u32) -> Vec<(u32, i64)> {
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
fn literal_rust_port_matches_c_extraction_term_for_term() {
    let cases = [
        (3_u32, 3_u32),
        (5, 5),
        (8, 8),
        (12, 12),
        (15, 15),
        (30, 30),
        (97, 1),
        (97, 97),
        (120, 1),
        (120, 120),
        (1009, 11),
        (2520, 26),
    ];
    let mut c = CContext::new().unwrap();
    let mut rust = RustContext::new();

    for (order, count) in cases {
        let lhs_terms = terms(order, count, 1);
        let rhs_terms = terms(order, count, 2);
        let c_lhs = c.from_terms(order, &lhs_terms).unwrap();
        let c_rhs = c.from_terms(order, &rhs_terms).unwrap();
        let rust_lhs = rust.from_terms(order, &lhs_terms).unwrap();
        let rust_rhs = rust.from_terms(order, &rhs_terms).unwrap();

        assert_eq!(rust_lhs.order(), c_lhs.order(), "lhs order {order}");
        assert_eq!(rust_lhs.terms(), c_lhs.terms(), "lhs terms {order}");
        assert_eq!(rust_rhs.order(), c_rhs.order(), "rhs order {order}");
        assert_eq!(rust_rhs.terms(), c_rhs.terms(), "rhs terms {order}");

        let c_sum = c.add(&c_lhs, &c_rhs).unwrap();
        let rust_sum = rust.add(&rust_lhs, &rust_rhs).unwrap();
        assert_eq!(rust_sum.order(), c_sum.order(), "sum order {order}");
        assert_eq!(rust_sum.terms(), c_sum.terms(), "sum terms {order}");

        let c_product = c.mul(&c_lhs, &c_rhs).unwrap();
        let rust_product = rust.mul(&rust_lhs, &rust_rhs).unwrap();
        assert_eq!(
            rust_product.order(),
            c_product.order(),
            "product order {order}"
        );
        assert_eq!(
            rust_product.terms(),
            c_product.terms(),
            "product terms {order}"
        );
    }
}

#[test]
fn literal_rust_port_matches_c_extraction_for_mixed_fields() {
    let pairs = [(3_u32, 4_u32), (3, 5), (4, 6), (5, 8), (8, 9), (12, 15)];
    let mut c = CContext::new().unwrap();
    let mut rust = RustContext::new();

    for (left_order, right_order) in pairs {
        let left_terms = terms(left_order, left_order.min(7), 3);
        let right_terms = terms(right_order, right_order.min(7), 4);
        let c_left = c.from_terms(left_order, &left_terms).unwrap();
        let c_right = c.from_terms(right_order, &right_terms).unwrap();
        let rust_left = rust.from_terms(left_order, &left_terms).unwrap();
        let rust_right = rust.from_terms(right_order, &right_terms).unwrap();

        let c_sum = c.add(&c_left, &c_right).unwrap();
        let rust_sum = rust.add(&rust_left, &rust_right).unwrap();
        assert_eq!(rust_sum.order(), c_sum.order());
        assert_eq!(rust_sum.terms(), c_sum.terms());

        let c_product = c.mul(&c_left, &c_right).unwrap();
        let rust_product = rust.mul(&rust_left, &rust_right).unwrap();
        assert_eq!(rust_product.order(), c_product.order());
        assert_eq!(rust_product.terms(), c_product.terms());
    }
}

#[cfg(feature = "libgap")]
#[test]
fn literal_rust_products_equal_unmodified_gap() {
    use gap_cyclotom_reference::libgap::Context as LibgapContext;

    let gap = LibgapContext::new().unwrap();
    let mut rust = RustContext::new();
    for (order, count) in [(5_u32, 5_u32), (12, 12), (97, 1), (120, 120), (1009, 11)] {
        let left_terms = terms(order, count, 1);
        let right_terms = terms(order, count, 2);
        let gap_left = gap.from_integer_terms(order, &left_terms).unwrap();
        let gap_right = gap.from_integer_terms(order, &right_terms).unwrap();
        let gap_product = gap.mul(&gap_left, &gap_right).unwrap();

        let rust_left = rust.from_terms(order, &left_terms).unwrap();
        let rust_right = rust.from_terms(order, &right_terms).unwrap();
        let rust_product = rust.mul(&rust_left, &rust_right).unwrap();
        let rust_as_gap = gap
            .from_integer_terms(rust_product.order(), &rust_product.terms())
            .unwrap();
        assert!(
            rust_as_gap == gap_product,
            "libgap mismatch for order {order}"
        );
    }
}
