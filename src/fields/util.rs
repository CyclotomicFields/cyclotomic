#[derive(Eq, PartialEq)]
pub enum Sign {
    Plus,
    Minus,
}

pub fn are_coprime(x: i64, y: i64) -> bool {
    num::integer::gcd(x as u64, y as u64) == 1
}

pub fn phi(n: i64) -> i64 {
    let mut count = 0;
    for k in 1..n {
        if are_coprime(n, k) {
            count += 1;
        }
    }
    count
}

pub fn math_mod(x: &i64, n: &i64) -> i64 {
    let res = (x % n + n) % n;
    res
}

pub fn count_powers(n: &i64, n_divisors: &Vec<i64>) -> Vec<(i64, i64)> {
    let mut result = vec![];
    let mut n_factored = n.clone();

    for divisor in n_divisors {
        let mut power: u64 = 0;

        while n_factored % divisor == 0 {
            power += 1;
            n_factored = n_factored / divisor;
        }

        if power != 0 {
            result.push((divisor.clone(), power as i64));
        }
    }

    result
}
