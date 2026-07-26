use gap_cyclotom_reference::libgap::{CharacterTableCase, Context as LibgapContext};
use gap_cyclotom_reference::tagged::Context as TaggedContext;
use gap_cyclotom_reference::tagged_character::PreparedCharacterTable;
use std::hint::black_box;
use std::time::{Duration, Instant};

struct Measurement {
    iterations: u64,
    median_ns: f64,
    minimum_ns: f64,
    maximum_ns: f64,
}

fn environment_u64(name: &str, default: u64) -> u64 {
    std::env::var(name)
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(default)
}

fn measure(mut operation: impl FnMut(), sample_time: Duration, samples: usize) -> Measurement {
    let mut iterations = 1_u64;
    loop {
        let start = Instant::now();
        for _ in 0..iterations {
            operation();
        }
        if start.elapsed() >= sample_time || iterations >= 1 << 20 {
            break;
        }
        iterations *= 2;
    }

    let mut observations = Vec::with_capacity(samples);
    for _ in 0..samples {
        let start = Instant::now();
        for _ in 0..iterations {
            operation();
        }
        observations.push(start.elapsed().as_nanos() as f64 / iterations as f64);
    }
    observations.sort_by(f64::total_cmp);
    Measurement {
        iterations,
        median_ns: observations[observations.len() / 2],
        minimum_ns: observations[0],
        maximum_ns: observations[observations.len() - 1],
    }
}

fn selected(value: &str, filter: &str) -> bool {
    filter == "all" || value.eq_ignore_ascii_case(filter)
}

fn main() {
    let sample_time = Duration::from_millis(environment_u64("CYCLOTOMIC_BENCH_SAMPLE_MS", 25));
    let samples = environment_u64("CYCLOTOMIC_BENCH_SAMPLES", 15) as usize;
    assert!(samples > 0, "CYCLOTOMIC_BENCH_SAMPLES must be positive");
    let table_filter = std::env::var("CYCLOTOMIC_BENCH_TABLE").unwrap_or_else(|_| "all".into());
    let implementation_filter =
        std::env::var("CYCLOTOMIC_BENCH_IMPL").unwrap_or_else(|_| "all".into());

    let gap = LibgapContext::new().expect("initialize unmodified libgap");
    let mut tagged = TaggedContext::new();
    println!(
        "character tensor decompositions: {samples} samples, target {} ms/sample",
        sample_time.as_millis()
    );

    for case in CharacterTableCase::ALL {
        if !selected(case.name(), &table_filter) {
            continue;
        }
        let table = gap
            .character_table(case)
            .expect("construct GAP character table");
        let prepared =
            PreparedCharacterTable::import(&mut tagged, &table).expect("import character table");
        let pairs: Vec<_> = (0..table.rows.len())
            .flat_map(|lhs| (lhs..table.rows.len()).map(move |rhs| (lhs, rhs)))
            .collect();

        for &(lhs, rhs) in &pairs {
            let expected = gap
                .character_tensor_decomposition(case, lhs, rhs, table.rows.len())
                .expect("decompose tensor product in GAP");
            let actual = prepared
                .tensor_multiplicities(&mut tagged, lhs, rhs)
                .expect("decompose tensor product in tagged Rust");
            assert_eq!(actual, expected, "{} rows {lhs}, {rhs}", case.name());
        }

        if selected("gap", &implementation_filter) {
            let measurement = measure(
                || {
                    for &(lhs, rhs) in &pairs {
                        black_box(
                            gap.character_tensor_decomposition(case, lhs, rhs, table.rows.len())
                                .unwrap(),
                        );
                    }
                },
                sample_time,
                samples,
            );
            print_measurement(case, "gap", pairs.len(), &measurement);
        }

        if selected("rust", &implementation_filter) {
            let measurement = measure(
                || {
                    for &(lhs, rhs) in &pairs {
                        black_box(
                            prepared
                                .tensor_multiplicities(&mut tagged, lhs, rhs)
                                .unwrap(),
                        );
                    }
                },
                sample_time,
                samples,
            );
            print_measurement(case, "rust", pairs.len(), &measurement);
        }
    }
}

fn print_measurement(
    case: CharacterTableCase,
    implementation: &str,
    pair_count: usize,
    measurement: &Measurement,
) {
    println!(
        "{:<10} {:<4} {:>9.3} us/tensor  median; range {:>9.3}..{:>9.3}; {} pairs × {} iterations/sample",
        case.name(),
        implementation,
        measurement.median_ns / pair_count as f64 / 1_000.0,
        measurement.minimum_ns / pair_count as f64 / 1_000.0,
        measurement.maximum_ns / pair_count as f64 / 1_000.0,
        pair_count,
        measurement.iterations,
    );
}
