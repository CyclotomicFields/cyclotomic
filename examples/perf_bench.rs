use cyclotomic::fields::dense;
use cyclotomic::fields::exponent::Exponent;
use cyclotomic::fields::sparse::{self, ExpCoeffMap};
use cyclotomic::fields::structure::CyclotomicField;
use cyclotomic::fields::{
    AdditiveGroupElement, CyclotomicFieldElement, MultiplicativeGroupElement, Q,
};
use std::env;
use std::hint::black_box;
use std::time::{Duration, Instant};

#[derive(Clone)]
struct Record {
    representation: &'static str,
    operation: &'static str,
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

fn dense_element(order: i64, density: f64, salt: i64) -> dense::Number {
    let mut coeffs = vec![Q::from(0); order as usize];
    for exp in 0..order {
        if include_term(exp + salt, order, density) {
            coeffs[exp as usize] = rational_for(order + salt, exp);
        }
    }
    dense::Number::new(&order, &coeffs)
}

fn sparse_element(order: i64, density: f64, salt: i64) -> sparse::Number<i64, Q> {
    let mut coeffs = ExpCoeffMap::default();
    for exp in 0..order {
        if include_term(exp + salt, order, density) {
            coeffs.insert(exp, rational_for(order + salt, exp));
        }
    }
    sparse::Number::<i64, Q>::new(&order, &coeffs)
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
    operation: &'static str,
    order: i64,
    density: f64,
    iterations: u64,
    ns_per_iter: f64,
) {
    records.push(Record {
        representation,
        operation,
        order,
        phi: Exponent::phi(&order),
        density,
        iterations,
        ns_per_iter,
        ops_per_second: 1_000_000_000.0 / ns_per_iter,
    });
}

fn benchmark_dense(
    records: &mut Vec<Record>,
    orders: &[i64],
    densities: &[f64],
    min_duration: Duration,
) {
    for &order in orders {
        eprintln!("benchmarking dense order {}", order);
        for &density in densities {
            let left = dense_element(order, density, 1);
            let right = dense_element(order, density, 2);
            let scalar = Q::from((7, 5));

            let (iterations, ns) = time_loop(min_duration, 1 << 24, || {
                let mut z1 = left.clone();
                let mut z2 = right.clone();
                black_box(z1.add(&mut z2));
            });
            push_record(records, "dense", "add", order, density, iterations, ns);

            let (iterations, ns) = time_loop(min_duration, 1 << 24, || {
                let mut z = left.clone();
                black_box(z.scalar_mul(&scalar));
            });
            push_record(
                records,
                "dense",
                "scalar_mul",
                order,
                density,
                iterations,
                ns,
            );

            let (iterations, ns) = time_loop(min_duration, 1 << 18, || {
                let mut z1 = left.clone();
                let mut z2 = right.clone();
                black_box(z1.mul(&mut z2));
            });
            push_record(records, "dense", "mul", order, density, iterations, ns);
        }
    }
}

fn benchmark_sparse(
    records: &mut Vec<Record>,
    orders: &[i64],
    densities: &[f64],
    min_duration: Duration,
) {
    for &order in orders {
        eprintln!("benchmarking sparse order {}", order);
        for &density in densities {
            let left = sparse_element(order, density, 1);
            let right = sparse_element(order, density, 2);
            let scalar = Q::from((7, 5));

            let (iterations, ns) = time_loop(min_duration, 1 << 24, || {
                let mut z1 = left.clone();
                let mut z2 = right.clone();
                black_box(z1.add(&mut z2));
            });
            push_record(records, "sparse", "add", order, density, iterations, ns);

            let (iterations, ns) = time_loop(min_duration, 1 << 24, || {
                let mut z = left.clone();
                black_box(z.scalar_mul(&scalar));
            });
            push_record(
                records,
                "sparse",
                "scalar_mul",
                order,
                density,
                iterations,
                ns,
            );

            let (iterations, ns) = time_loop(min_duration, 1 << 18, || {
                let mut z1 = left.clone();
                let mut z2 = right.clone();
                black_box(z1.mul(&mut z2));
            });
            push_record(records, "sparse", "mul", order, density, iterations, ns);
        }
    }
}

fn benchmark_structure(
    records: &mut Vec<Record>,
    orders: &[i64],
    densities: &[f64],
    min_duration: Duration,
) {
    for &order in orders {
        eprintln!("benchmarking structure order {}", order);
        let (iterations, ns) = time_loop(min_duration, 1 << 14, || {
            black_box(CyclotomicField::new(order));
        });
        push_record(
            records,
            "structure",
            "construct",
            order,
            1.0,
            iterations,
            ns,
        );

        let field = CyclotomicField::new(order);
        for &density in densities {
            let left = structure_element(&field, order, density, 1);
            let right = structure_element(&field, order, density, 2);

            let (iterations, ns) = time_loop(min_duration, 1 << 18, || {
                black_box(field.mul(&left, &right));
            });
            push_record(records, "structure", "mul", order, density, iterations, ns);

            let (iterations, ns) = time_loop(min_duration, 1 << 18, || {
                black_box(field.mul_flat(&left, &right));
            });
            push_record(
                records,
                "structure_flat",
                "mul",
                order,
                density,
                iterations,
                ns,
            );

            let (iterations, ns) = time_loop(min_duration, 1 << 18, || {
                black_box(field.mul_sparse_constants(&left, &right));
            });
            push_record(
                records,
                "structure_sparse",
                "mul",
                order,
                density,
                iterations,
                ns,
            );

            let (iterations, ns) = time_loop(min_duration, 1 << 18, || {
                black_box(field.mul_sparse_constants_skip_zero_inputs(&left, &right));
            });
            push_record(
                records,
                "structure_sparse_skip_zero",
                "mul",
                order,
                density,
                iterations,
                ns,
            );
        }
    }
}

fn print_json(records: &[Record]) {
    println!("[");
    for (i, record) in records.iter().enumerate() {
        let comma = if i + 1 == records.len() { "" } else { "," };
        println!(
            "  {{\"representation\":\"{}\",\"operation\":\"{}\",\"order\":{},\"phi\":{},\"density\":{:.2},\"iterations\":{},\"ns_per_iter\":{:.3},\"ops_per_second\":{:.3}}}{}",
            record.representation,
            record.operation,
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
    let quick = env::args().any(|arg| arg == "--quick");
    let min_duration = if quick {
        Duration::from_millis(20)
    } else {
        Duration::from_millis(175)
    };
    let orders = if quick {
        vec![5, 8, 12, 15]
    } else {
        vec![
            // Small primes.
            5, 7, 11, 13, 17, 19, // Larger primes up to 100.
            23, 29, 31, 37, 43, 53, 61, 71, 83, 97, // Powers of two.
            4, 8, 16, 32, 64, // Highly factorizable and mixed composite orders.
            12, 18, 24, 30, 36, 40, 48, 60, 72, 84, 90, 96, 100,
        ]
    };
    let structure_orders = if quick {
        vec![5, 8, 12]
    } else {
        vec![5, 7, 8, 12, 16, 24, 30, 32, 40]
    };
    let densities = vec![0.10, 0.25, 0.50, 0.75, 1.00];

    let mut records = Vec::new();
    benchmark_dense(&mut records, &orders, &densities, min_duration);
    benchmark_sparse(&mut records, &orders, &densities, min_duration);
    benchmark_structure(&mut records, &structure_orders, &densities, min_duration);
    print_json(&records);
}
