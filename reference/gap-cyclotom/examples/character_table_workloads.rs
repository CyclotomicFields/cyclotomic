use cyclotomic::fields::dense;
use cyclotomic::fields::rational::{HybridRational, Rational as RationalCoefficient};
use cyclotomic::fields::sparse;
use cyclotomic::fields::{AdditiveGroupElement, CyclotomicFieldElement, GenericCyclotomic};
use gap_cyclotom_reference::libgap::{CharacterTable, CharacterTableCase, Context, Cyclotomic};
use rug::Integer;
use std::hint::black_box;
use std::time::{Duration, Instant};

fn generic(value: &Cyclotomic) -> GenericCyclotomic {
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

fn convert_table<T, Q>(table: &CharacterTable) -> Vec<Vec<T>>
where
    T: CyclotomicFieldElement<i64, Q>,
    Q: RationalCoefficient,
{
    table
        .rows
        .iter()
        .map(|row| {
            row.iter()
                .map(|value| T::from_generic(&generic(value)))
                .collect()
        })
        .collect()
}

fn tensor_decomposition<T, Q>(
    rows: &[Vec<T>],
    weighted_conjugates: &[Vec<T>],
    group_order: i64,
    lhs: usize,
    rhs: usize,
) -> Vec<T>
where
    T: CyclotomicFieldElement<i64, Q>,
    Q: RationalCoefficient,
{
    let mut sums: Vec<Option<T>> = vec![None; weighted_conjugates.len()];
    for class in 0..rows[lhs].len() {
        let mut product = rows[lhs][class].clone();
        product.mul(&mut rows[rhs][class].clone());
        for (sum, irreducible) in sums.iter_mut().zip(weighted_conjugates) {
            let mut term = product.clone();
            term.mul(&mut irreducible[class].clone());
            if let Some(sum) = sum {
                sum.add(&mut term);
            } else {
                *sum = Some(term);
            }
        }
    }
    sums.into_iter()
        .map(|sum| {
            let mut sum = sum.expect("character tables have at least one class");
            sum.scalar_mul(&Q::from((1, group_order as u64)));
            sum
        })
        .collect()
}

fn weighted_conjugates<T, Q>(rows: &[Vec<T>], class_sizes: &[i64]) -> Vec<Vec<T>>
where
    T: CyclotomicFieldElement<i64, Q>,
    Q: RationalCoefficient,
{
    rows.iter()
        .map(|row| {
            row.iter()
                .zip(class_sizes)
                .map(|(character, class_size)| {
                    let mut weighted = character.complex_conjugate();
                    weighted.scalar_mul(&Q::from(*class_size));
                    weighted
                })
                .collect()
        })
        .collect()
}

fn tensor_decomposition_packed(
    rows: &[Vec<sparse::Number<i64, HybridRational>>],
    weighted_conjugates: &[Vec<sparse::Number<i64, HybridRational>>],
    group_order: i64,
    lhs: usize,
    rhs: usize,
    scratch: &mut sparse::mul::PackedMulScratch<HybridRational>,
) -> Vec<sparse::Number<i64, HybridRational>> {
    let mut sums: Vec<Option<sparse::Number<i64, HybridRational>>> =
        vec![None; weighted_conjugates.len()];
    for class in 0..rows[lhs].len() {
        let product = rows[lhs][class].mul_packed(&rows[rhs][class], scratch);
        for (sum, irreducible) in sums.iter_mut().zip(weighted_conjugates) {
            let mut term = product.mul_packed(&irreducible[class], scratch);
            if let Some(sum) = sum {
                sum.add(&mut term);
            } else {
                *sum = Some(term);
            }
        }
    }
    sums.into_iter()
        .map(|sum| {
            let mut sum = sum.expect("character tables have at least one class");
            sum.scalar_mul(&HybridRational::from((1, group_order as u64)));
            sum
        })
        .collect()
}

fn sparse_multiplicities(
    rows: &[Vec<sparse::Number<i64, HybridRational>>],
    weighted_conjugates: &[Vec<sparse::Number<i64, HybridRational>>],
    group_order: i64,
    lhs: usize,
    rhs: usize,
    scratch: &mut sparse::mul::PackedMulScratch<HybridRational>,
) -> Vec<HybridRational> {
    tensor_decomposition_packed(
        rows,
        weighted_conjugates,
        group_order,
        lhs,
        rhs,
        scratch,
    )
    .into_iter()
    .map(|value| {
        value
            .into_rational()
            .expect("character inner products are rational")
    })
    .collect()
}

fn dense_multiplicities(
    rows: &[Vec<dense::Number>],
    weighted_conjugates: &[Vec<dense::Number>],
    group_order: i64,
    lhs: usize,
    rhs: usize,
) -> Vec<rug::Rational> {
    tensor_decomposition(
        rows,
        weighted_conjugates,
        group_order,
        lhs,
        rhs,
    )
    .into_iter()
    .map(|value| {
        value
            .into_rational()
            .expect("character inner products are rational")
    })
    .collect()
}

fn assert_rational_decomposition<Q: RationalCoefficient>(actual: &[Q], expected: &[i64]) {
    assert_eq!(actual.len(), expected.len());
    for (actual, expected) in actual.iter().zip(expected) {
        assert_eq!(actual.numer(), Integer::from(*expected));
        assert_eq!(actual.denom(), Integer::from(1));
    }
}

fn interesting_rows(table: &CharacterTable) -> (usize, usize) {
    let rows: Vec<_> = table
        .rows
        .iter()
        .enumerate()
        .filter(|(_, row)| {
            row.iter()
                .any(|value| value.order().expect("read GAP conductor") > 1)
        })
        .map(|(index, _)| index)
        .collect();
    match rows.as_slice() {
        [first, second, ..] => (*first, *second),
        [first] => (*first, *first),
        [] => {
            let last = table.rows.len() - 1;
            (last.saturating_sub(1), last)
        }
    }
}

fn measure(mut operation: impl FnMut()) -> (u64, f64) {
    let minimum = Duration::from_millis(10);
    let mut iterations = 1_u64;
    loop {
        let start = Instant::now();
        for _ in 0..iterations {
            operation();
        }
        let elapsed = start.elapsed();
        if elapsed >= minimum || iterations >= 1 << 14 {
            return (iterations, elapsed.as_nanos() as f64 / iterations as f64);
        }
        iterations *= 2;
    }
}

fn main() {
    let context = Context::new().expect("initialize unmodified libgap");
    println!("[");
    let mut first_record = true;

    for case in CharacterTableCase::ALL {
        let table = context
            .character_table(case)
            .expect("construct GAP character table");
        let (lhs, rhs) = interesting_rows(&table);
        let expected = context
            .character_tensor_decomposition(case, lhs, rhs, table.rows.len())
            .expect("decompose tensor product in GAP");
        let sparse_rows: Vec<Vec<sparse::Number<i64, HybridRational>>> =
            convert_table(&table);
        let dense_rows: Vec<Vec<dense::Number>> = convert_table(&table);
        let sparse_weighted_conjugates =
            weighted_conjugates(&sparse_rows, &table.class_sizes);
        let dense_weighted_conjugates = weighted_conjugates(&dense_rows, &table.class_sizes);
        let mut sparse_scratch = sparse::mul::PackedMulScratch::new();

        assert_rational_decomposition(
            &sparse_multiplicities(
                &sparse_rows,
                &sparse_weighted_conjugates,
                table.group_order,
                lhs,
                rhs,
                &mut sparse_scratch,
            ),
            &expected,
        );
        assert_rational_decomposition(
            &dense_multiplicities(
                &dense_rows,
                &dense_weighted_conjugates,
                table.group_order,
                lhs,
                rhs,
            ),
            &expected,
        );

        let (gap_iterations, gap_ns) = measure(|| {
            black_box(
                context
                    .character_tensor_decomposition(case, lhs, rhs, table.rows.len())
                    .unwrap(),
            );
        });
        let (sparse_iterations, sparse_ns) = measure(|| {
            black_box(sparse_multiplicities(
                &sparse_rows,
                &sparse_weighted_conjugates,
                table.group_order,
                lhs,
                rhs,
                &mut sparse_scratch,
            ));
        });
        let (dense_iterations, dense_ns) = measure(|| {
            black_box(dense_multiplicities(
                &dense_rows,
                &dense_weighted_conjugates,
                table.group_order,
                lhs,
                rhs,
            ));
        });
        let (sparse_alloc_iterations, sparse_alloc_ns) = measure(|| {
            let rows: Vec<Vec<sparse::Number<i64, HybridRational>>> =
                convert_table(&table);
            let conjugates = weighted_conjugates(&rows, &table.class_sizes);
            let mut scratch = sparse::mul::PackedMulScratch::new();
            black_box(sparse_multiplicities(
                &rows,
                &conjugates,
                table.group_order,
                lhs,
                rhs,
                &mut scratch,
            ));
        });
        let (dense_alloc_iterations, dense_alloc_ns) = measure(|| {
            let rows: Vec<Vec<dense::Number>> = convert_table(&table);
            let conjugates = weighted_conjugates(&rows, &table.class_sizes);
            black_box(dense_multiplicities(
                &rows,
                &conjugates,
                table.group_order,
                lhs,
                rhs,
            ));
        });

        for (implementation, mode, iterations, ns) in [
            (
                "gap_unmodified_libgap",
                "native_representation",
                gap_iterations,
                gap_ns,
            ),
            (
                "rust_sparse_packed_hybrid",
                "prepared_arithmetic_only",
                sparse_iterations,
                sparse_ns,
            ),
            (
                "rust_dense_rational",
                "prepared_arithmetic_only",
                dense_iterations,
                dense_ns,
            ),
            (
                "rust_sparse_packed_hybrid",
                "allocation_inclusive_end_to_end",
                sparse_alloc_iterations,
                sparse_alloc_ns,
            ),
            (
                "rust_dense_rational",
                "allocation_inclusive_end_to_end",
                dense_alloc_iterations,
                dense_alloc_ns,
            ),
        ] {
            if !first_record {
                println!(",");
            }
            first_record = false;
            print!(
                "  {{\"family\":\"character_table\",\"table\":\"{}\",\"implementation\":\"{implementation}\",\"operation\":\"tensor_decomposition\",\"mode\":\"{mode}\",\"rows\":{},\"classes\":{},\"group_order\":{},\"lhs_row\":{lhs},\"rhs_row\":{rhs},\"iterations\":{iterations},\"ns_per_iter\":{ns:.3}}}",
                case.name(),
                table.rows.len(),
                table.class_sizes.len(),
                table.group_order,
            );
        }
    }
    println!("\n]");
}
