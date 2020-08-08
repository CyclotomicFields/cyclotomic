use crate::fields::Z;

#[derive(Eq, PartialEq)]
pub enum Sign {
    Plus,
    Minus,
}

pub fn are_coprime(x: i64, y: i64) -> bool {
    num::integer::gcd(x as u64, y as u64) == 1
}

// TODO: this is terrible, use one of rob's good ones
pub fn phi(n: i64) -> i64 {
    let mut count = 0;
    for k in 1..n {
        if are_coprime(n, k) {
            count += 1;
        }
    }
    count
}

// TODO: make generic, combine with phi
pub fn phi_big(n: &Z) -> Z {
    let mut count = Z::from(0);
    let mut k = Z::from(1);
    while &k != n {
        if n.gcd_ref(&k).into(): Z == 1 {
            count += 1;
        }
        k += 1;
    }
    count
}

// TODO: make generic, combine with math_mod
pub fn math_mod_big(x: &Z, n: &Z) -> Z {
    let rem: Z = (x % n).into();
    let add: Z = (rem + n).into();
    add % n
}

pub fn math_mod(x: &i64, n: &i64) -> i64 {
    (x % n + n) % n
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

pub fn count_powers_big(n: &Z, n_divisors: &Vec<Z>) -> Vec<(Z, u32)> {
    let mut result = vec![];
    let mut n_factored = n.clone();

    for divisor in n_divisors {
        let mut power = 0;

        while (&n_factored % divisor).into(): Z == 0 {
            power += 1;
            n_factored = (&n_factored / divisor).into();
        }

        if power != 0 {
            result.push((divisor.clone(), power));
        }
    }

    result
}
