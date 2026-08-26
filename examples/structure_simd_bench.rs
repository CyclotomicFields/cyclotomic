use cyclotomic::fields::exponent::Exponent;
use cyclotomic::fields::structure::CyclotomicField;
use cyclotomic::fields::Q;
use std::hint::black_box;
use std::time::{Duration, Instant};

struct Record {
    representation: &'static str,
    order: i64,
    phi: i64,
    density: f64,
    iterations: u64,
    ns_per_iter: f64,
    ops_per_second: f64,
}

fn rational_for(order: i64, index: i64) -> Q {
    let numerator = ((order * 17 + index * 31) % 23) - 11;
    let denominator = (((order + index * 7).abs() % 9) + 1) as u64;
    Q::from((numerator, denominator))
}

fn include_term(index: i64, order: i64, density: f64) -> bool {
    if density >= 1.0 {
        return true;
    }
    let bucket = ((index * 37 + order * 11) % 100) as f64 / 100.0;
    bucket < density
}

fn structure_element(field: &CyclotomicField, order: i64, density: f64, salt: i64) -> Vec<Q> {
    let mut result = vec![Q::from(0); field.basis.len()];
    for (i, coeff) in result.iter_mut().enumerate() {
        if include_term(i as i64 + salt, order, density) {
            *coeff = rational_for(order + salt, i as i64);
        }
    }
    result
}

fn time_loop<F>(min_duration: Duration, max_iters: u64, mut f: F) -> (u64, f64)
where
    F: FnMut(),
{
    let mut iterations = 1_u64;
    loop {
        let start = Instant::now();
        for _ in 0..iterations {
            f();
        }
        let elapsed = start.elapsed();
        if elapsed >= min_duration || iterations >= max_iters {
            return (iterations, elapsed.as_nanos() as f64 / iterations as f64);
        }
        iterations = (iterations * 2).min(max_iters);
    }
}

fn push_record(
    records: &mut Vec<Record>,
    representation: &'static str,
    order: i64,
    density: f64,
    iterations: u64,
    ns_per_iter: f64,
) {
    records.push(Record {
        representation,
        order,
        phi: Exponent::phi(&order),
        density,
        iterations,
        ns_per_iter,
        ops_per_second: 1_000_000_000.0 / ns_per_iter,
    });
}

fn print_json(records: &[Record]) {
    println!("[");
    for (i, record) in records.iter().enumerate() {
        let comma = if i + 1 == records.len() { "" } else { "," };
        println!(
            "  {{\"representation\":\"{}\",\"operation\":\"mul\",\"order\":{},\"phi\":{},\"density\":{:.2},\"iterations\":{},\"ns_per_iter\":{:.3},\"ops_per_second\":{:.3}}}{}",
            record.representation,
            record.order,
            record.phi,
            record.density,
            record.iterations,
            record.ns_per_iter,
            record.ops_per_second,
            comma
        );
    }
    println!("]");
}

fn main() {
    let min_duration = Duration::from_millis(150);
    let orders = [8, 12, 16, 24, 30, 32, 40];
    let densities = [0.10, 0.50, 1.00];
    let mut records = Vec::new();

    for order in orders {
        eprintln!("benchmarking structure order {}", order);
        let field = CyclotomicField::new(order);
        for density in densities {
            let left = structure_element(&field, order, density, 1);
            let right = structure_element(&field, order, density, 2);
            let left_i64 = field.to_i64_shared_den_element(&left).unwrap();
            let right_i64 = field.to_i64_shared_den_element(&right).unwrap();

            let (iterations, ns) = time_loop(min_duration, 1 << 18, || {
                black_box(field.mul(&left, &right));
            });
            push_record(&mut records, "structure", order, density, iterations, ns);

            let (iterations, ns) = time_loop(min_duration, 1 << 18, || {
                black_box(field.mul_sparse_constants_skip_zero_inputs(&left, &right));
            });
            push_record(
                &mut records,
                "structure_sparse_skip_zero",
                order,
                density,
                iterations,
                ns,
            );

            let (iterations, ns) = time_loop(min_duration, 1 << 18, || {
                black_box(field.mul_i64_flat_shared_den(&left, &right).unwrap());
            });
            push_record(
                &mut records,
                "structure_i64_flat_shared_den",
                order,
                density,
                iterations,
                ns,
            );

            let (iterations, ns) = time_loop(min_duration, 1 << 18, || {
                black_box(field.mul_i64_sparse_shared_den(&left, &right).unwrap());
            });
            push_record(
                &mut records,
                "structure_i64_sparse_shared_den",
                order,
                density,
                iterations,
                ns,
            );

            let (iterations, ns) = time_loop(min_duration, 1 << 18, || {
                black_box(field.mul_i64_output_grouped_shared_den(&left, &right).unwrap());
            });
            push_record(
                &mut records,
                "structure_i64_output_grouped_shared_den",
                order,
                density,
                iterations,
                ns,
            );

            let (iterations, ns) = time_loop(min_duration, 1 << 18, || {
                black_box(
                    field
                        .mul_i64_sparse_preconverted(&left_i64, &right_i64)
                        .unwrap(),
                );
            });
            push_record(
                &mut records,
                "structure_i64_sparse_preconverted",
                order,
                density,
                iterations,
                ns,
            );

            let (iterations, ns) = time_loop(min_duration, 1 << 18, || {
                black_box(
                    field
                        .mul_i64_output_grouped_preconverted(&left_i64, &right_i64)
                        .unwrap(),
                );
            });
            push_record(
                &mut records,
                "structure_i64_output_grouped_preconverted",
                order,
                density,
                iterations,
                ns,
            );
        }
    }

    print_json(&records);
}
