use cyclotomic::fields::dense;
use cyclotomic::fields::rational::{HybridRational, Rational as RationalCoefficient};
use cyclotomic::fields::sparse;
use cyclotomic::fields::{CyclotomicFieldElement, GenericCyclotomic};
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
    class_sizes: &[i64],
    group_order: i64,
    lhs: usize,
    rhs: usize,
) -> Vec<T>
where
    T: CyclotomicFieldElement<i64, Q>,
    Q: RationalCoefficient,
{
    let product: Vec<_> = rows[lhs]
        .iter()
        .zip(&rows[rhs])
        .map(|(left, right)| {
            let mut result = left.clone();
            result.mul(&mut right.clone());
            result
        })
        .collect();

    rows.iter()
        .map(|irreducible| {
            let mut sum: Option<T> = None;
            for ((value, character), class_size) in product.iter().zip(irreducible).zip(class_sizes)
            {
                let mut term = value.clone();
                term.mul(&mut character.complex_conjugate());
                term.scalar_mul(&Q::from(*class_size));
                if let Some(sum) = &mut sum {
                    sum.add(&mut term);
                } else {
                    sum = Some(term);
                }
            }
            let mut sum = sum.expect("character tables have at least one class");
            sum.scalar_mul(&Q::from((1, group_order as u64)));
            sum
        })
        .collect()
}

fn integer(value: i64) -> GenericCyclotomic {
    GenericCyclotomic {
        order: Integer::from(1),
        exp_coeffs: [(Integer::from(0), (value, 1))].into_iter().collect(),
    }
}

fn assert_decomposition<T, Q>(actual: &[T], expected: &[i64])
where
    T: CyclotomicFieldElement<i64, Q>,
    Q: RationalCoefficient,
{
    assert_eq!(actual.len(), expected.len());
    for (actual, expected) in actual.iter().zip(expected) {
        let mut actual = actual.clone();
        let mut expected = T::from_generic(&integer(*expected));
        assert!(actual.eq(&mut expected));
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

        assert_decomposition(
            &tensor_decomposition(
                &sparse_rows,
                &table.class_sizes,
                table.group_order,
                lhs,
                rhs,
            ),
            &expected,
        );
        assert_decomposition(
            &tensor_decomposition(&dense_rows, &table.class_sizes, table.group_order, lhs, rhs),
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
            black_box(tensor_decomposition(
                &sparse_rows,
                &table.class_sizes,
                table.group_order,
                lhs,
                rhs,
            ));
        });
        let (dense_iterations, dense_ns) = measure(|| {
            black_box(tensor_decomposition(
                &dense_rows,
                &table.class_sizes,
                table.group_order,
                lhs,
                rhs,
            ));
        });

        for (implementation, iterations, ns) in [
            ("gap_unmodified_libgap", gap_iterations, gap_ns),
            ("rust_sparse_hybrid_rational", sparse_iterations, sparse_ns),
            ("rust_dense_rational", dense_iterations, dense_ns),
        ] {
            if !first_record {
                println!(",");
            }
            first_record = false;
            print!(
                "  {{\"family\":\"character_table\",\"table\":\"{}\",\"implementation\":\"{implementation}\",\"operation\":\"tensor_decomposition\",\"mode\":\"native_representation\",\"rows\":{},\"classes\":{},\"group_order\":{},\"lhs_row\":{lhs},\"rhs_row\":{rhs},\"iterations\":{iterations},\"ns_per_iter\":{ns:.3}}}",
                case.name(),
                table.rows.len(),
                table.class_sizes.len(),
                table.group_order,
            );
        }
    }
    println!("\n]");
}
