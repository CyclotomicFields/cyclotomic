use gap_cyclotom_reference::libgap::{CharacterTableCase, Context as LibgapContext};
use gap_cyclotom_reference::tagged::Context as TaggedContext;
use gap_cyclotom_reference::tagged_character::PreparedCharacterTable;
use std::hint::black_box;
use std::time::{Duration, Instant};

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
    let gap = LibgapContext::new().expect("initialize unmodified libgap");
    let mut tagged = TaggedContext::new();
    println!("[");
    let mut first = true;

    for case in CharacterTableCase::ALL {
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
            assert_eq!(actual, expected);
        }

        let (gap_iterations, gap_ns) = measure(|| {
            for &(lhs, rhs) in &pairs {
                black_box(
                    gap.character_tensor_decomposition(case, lhs, rhs, table.rows.len())
                        .unwrap(),
                );
            }
        });
        let (tagged_iterations, tagged_ns) = measure(|| {
            for &(lhs, rhs) in &pairs {
                black_box(
                    prepared
                        .tensor_multiplicities(&mut tagged, lhs, rhs)
                        .unwrap(),
                );
            }
        });

        for (implementation, mode, iterations, ns) in [
            (
                "gap_unmodified_libgap",
                "native_character_table",
                gap_iterations,
                gap_ns,
            ),
            (
                "rust_gap_port_tagged",
                "prepared_weighted_conjugates",
                tagged_iterations,
                tagged_ns,
            ),
        ] {
            if !first {
                println!(",");
            }
            first = false;
            print!(
                "  {{\"family\":\"character_table\",\"table\":\"{}\",\"implementation\":\"{implementation}\",\"operation\":\"all_unordered_tensor_decompositions\",\"mode\":\"{mode}\",\"rows\":{},\"classes\":{},\"tensor_pairs\":{},\"iterations\":{iterations},\"ns_per_sweep\":{ns:.3},\"ns_per_tensor\":{:.3}}}",
                case.name(),
                table.rows.len(),
                table.class_sizes.len(),
                pairs.len(),
                ns / pairs.len() as f64,
            );
        }
    }
    println!("\n]");
}
