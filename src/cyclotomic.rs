// \documentclass{article}
// \usepackage{amsfonts}
// \usepackage{amsmath}
// \usepackage{parskip}
// \usepackage[margin=1in]{geometry}
//
// \title{High-performance dense and sparse implementations of
// cyclotomic field arithmetic}
//
// \author{Kaashif Hymabaccus, Robert Moore}
// \date{Some time in 2020}
//
// \begin{document}
// \maketitle
//
// \begin{abstract}
//
// This paper presents a library for arithmetic using elements of the
// fields $\mathbb{Q}(\zeta_n)$, for $n < 2^{64}$. Existing
// implementations of this are usually part of a computer algebra
// system and thus not standalone. Examples are the implementations
// included with GAP and Magma. This makes it difficult to compare
// various encodings of elements of cyclotomic fields, and algorithms
// for the field operations. We implement the algorithms of GAP and
// Magma, describe several original algorithms, and present detailed
// benchmarks.
//
// \end{abstract}
//
// \section{Introduction}
//
// TODO: write something about GAP, Magma, cite some of those papers
//       about sparse and dense representations.
//
// \section{Implementation of GAP's algorithms}
//
// Although we are implementing this library in Rust, not C, we are
// aiming to match GAP's implementation as closely as possible, for
// comparison purposes with later encodings and algorithms.
//
// Suppose $z \in \mathbb{Q}(\zeta_n)$ where $\zeta_n$ is an $n$th
// root of unity. We can write $z = \sum_{i=0}^{n-1} c_i
// \zeta_n^{e_i}$ for some $c_i \in \mathbb{Q}$, $e_i \in \mathbb{Z}$.
//
// If we suppose $n$ is minimal, then GAP encodes $z$ as $(n, c, e)$
// where $n \in \mathbb{Z}_{>0}$, $c \in \mathbb{Q}^N$, $e \in
// [0..n-1]^N$, and $0 < N \leq n$. Additionally, the vector $e$ is
// sorted and duplicate-free. That is, if $i < j$, then $e_i <
// e_j$. This imposes a unique ordering on the terms of the sum. In
// essence, we have written $z$ as an element of the quotient ring
// $\mathbb{Q}[x]/(\Phi_n(x))$, while imposing an order on the terms
// of the polynomial.
extern crate num;

use std::cmp::Eq;
use std::ops::{Add, Mul};
use std::vec::Vec;
use self::num::{One, Zero};
use std::ptr::eq;

type Z = num::bigint::BigInt;
type Q = num::rational::BigRational;

#[derive(Debug, Clone)]
pub struct Cyclotomic {
    order: u64     /*   n   */,
    coeffs: Vec<Q> /*   c   */,
    exps: Vec<u64> /*   e   */,
}

impl Cyclotomic {
    pub fn new(order: u64, coeffs: Vec<Q>, exps: Vec<u64>) -> Cyclotomic {
        Cyclotomic {
            order,
            coeffs,
            exps,
        }
    }
}

// Without imposing further restrictions, the encoding of $z$ is not
// unique. Since we have chosen $n$ minimal, for uniqueness, we need
// only choose a basis for $K = \mathbb{Q}(\zeta_n)$. $K/\mathbb{Q}$
// is a degree $\phi(n)$ extension, and there is a natural choice of
// basis: $\{ \zeta_n^i : 0 \leq i \leq \phi(n)-1 \}$. This basis
// proves difficult to compute with. The choice of basis used in GAP
// is the set of $\zeta_n^i$ such that $i \notin
// (n/q)*[-(q/p-1)/2..(q/p-1)/2]$ mod $q$, for every odd prime divisor
// $p$ of $n$, where $q$ is the maximal power of $p$ in $n$, and $i
// \notin (n/q)*[q/2..q-1]$, if $q$ is the maximal power of 2 in $n$.
// It is not too difficult to see, that this gives in fact $\phi(n)$
// roots.

// TODO: I'm pretty sure it's very difficult to see. Maybe give a proof.

// TODO: tests!!! aaah

// TODO: heavily comment and explain this stuff, prove correctness etc

// TODO: implement reduction so the order doesn't blow up

impl Zero for Cyclotomic {
    fn zero() -> Self {
        Cyclotomic::new(1, vec![Q::new(Z::from(0), Z::from(1))], vec![0])
    }

    fn set_zero(&mut self) {
        self.order = 1;
        self.coeffs = vec![Q::new(Z::from(0), Z::from(1))];
        self.exps = vec![]
    }

    fn is_zero(&self) -> bool {
        return self.order == 1
            && self.coeffs == vec![Q::new(Z::from(0), Z::from(1))]
            && self.exps == vec![];
    }
}

impl One for Cyclotomic {
    fn one() -> Self {
        Cyclotomic::new(1, vec![Q::new(Z::from(1), Z::from(1))], vec![0])
    }
}

impl Cyclotomic {
    // TODO: right now we INCREASE the orders to get compatibility... we should try to
    // DECREASE order first!
    pub fn increase_order_to(z: &mut Cyclotomic, new_order: u64) -> () {
        z.exps = z
            .exps
            .clone()
            .into_iter()
            .map(|k| new_order * k / z.order)
            .collect();

        z.order = new_order
    }

    // TODO: Use Rob's code to actually do some reductions here?
    pub fn match_orders(z1: &mut Cyclotomic, z2: &mut Cyclotomic) -> () {
        let new_order = num::integer::lcm(z1.order, z2.order);
        Cyclotomic::increase_order_to(z1, new_order);
        Cyclotomic::increase_order_to(z2, new_order);
    }
}

impl PartialEq for Cyclotomic {
    // Note: checking equality mutates the encoding of $z_1$ and $z_2$, but does not
    // change the values represented. Since we mutate, we cannot use the Rust traits
    // \tt{Eq} and \tt{PartialEq}. The same is true for all arithmetic operations.
    fn eq(&self, other: &Self) -> bool {
        let mut z1 = self.clone();
        let mut z2 = other.clone();
        Cyclotomic::match_orders(&mut z1, &mut z2);
        // Since we assume the exponent list is sorted, we can just do a
        // comparison of the coeff and exp vectors, the number is in
        // canonical form.
        z1.coeffs == z2.coeffs && z1.exps == z2.exps
    }

    fn ne(&self, other: &Self) -> bool {
        return !eq(self, other);
    }
}

impl Add for Cyclotomic {
    type Output = Cyclotomic;

    fn add(self, rhs: Self) -> Self::Output {
        let mut z1 = self.clone();
        let mut z2 = rhs.clone();
        Cyclotomic::match_orders(&mut z1, &mut z2);

        let mut z1_i = 0;
        let mut z2_i = 0;
        let mut result = Cyclotomic::zero();

        // We are producing a result cyclotomic in canonical form assuming the
        // inputs are in canonical form. That is, we build the result term by
        // term, assuming that the exponent lists are sorted.
        while z1_i < z1.exps.len() || z2_i < z2.exps.len() {
            while z1_i < z1.exps.len() && z1.exps[z1_i] < z2.exps[z2_i] {
                result.coeffs.push(z1.coeffs[z1_i].clone());
                result.exps.push(z1.exps[z1_i].clone());
                z1_i += 1;
            }

            while z2_i < z2.exps.len() && z2.exps[z2_i] < z1.exps[z1_i] {
                result.coeffs.push(z2.coeffs[z2_i].clone());
                result.exps.push(z2.exps[z2_i].clone());
                z2_i += 1;
            }

            // If the exponents match up, we just add the coefficients - they are for the same term
            if z1.exps[z1_i] == z2.exps[z2_i] {
                result
                    .coeffs
                    .push(z1.coeffs[z1_i].clone() + z2.coeffs[z2_i].clone());
                result.exps.push(z1.exps[z1_i]);
                z1_i += 1;
                z2_i += 1
            }
        }

        result
    }
}

impl Mul for Cyclotomic {
    type Output = Cyclotomic;

    // I could have a guess at how to implement this, but I would inevitably embarrass myself.
    fn mul(self, rhs: Self) -> Self::Output {
        unimplemented!()
    }
}

// \section{Improving performance}
//
// Insert our incredible mathematical, algorithmic genius here.

// \end{document}
