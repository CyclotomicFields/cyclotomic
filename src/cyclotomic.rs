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
// about sparse and dense representations.
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
use std::ops::Add;
use std::vec::Vec;

type Z = num::bigint::BigInt;
type Q = num::rational::BigRational;

#[derive(Debug)]
struct Cyclotomic {
    order: u64, // n
    coeffs: Vec<Q>, // c
    exps: Vec<u64>, // e
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

// TODO: I'm pretty sure it's very difficult to see. Maybe give a
// proof.

// TODO: tests!!! aaah

// TODO: heavily comment and explain this stuff, prove correctness etc

// TODO: implement reduction so the order doesn't blow up

impl Cyclotomic {
    pub fn new(order: u64, coeffs: Vec<Q>, exps: Vec<u64>) -> Cyclotomic {
        Cyclotomic {
            order,
            coeffs,
            exps,
        }
    }
    pub fn zero() -> Cyclotomic {
        Cyclotomic::new(1, vec![Q::new(Z::from(0), Z::from(1))], vec![0])
    }
    pub fn one() -> Cyclotomic {
        Cyclotomic::new(1, vec![Q::new(Z::from(1), Z::from(1))], vec![0])
    }
    // Expresses self as a sum of powers of the Nth root of unity
    pub fn increase_order_to(&self, new_order: u64) -> Cyclotomic {
        Cyclotomic::new(
            new_order,
            self.coeffs.clone(),
            self.exps
                .clone()
                .into_iter()
                .map(|k| new_order * k / self.order)
                .collect(),
        )
    }

    // Returns a pair of cyclotomics representing the same numbers, but
    // expressed as elements of the same cyclotomic field
    pub fn match_orders(z1: &Cyclotomic, z2: &Cyclotomic) -> (Cyclotomic, Cyclotomic) {
        let new_order = num::integer::lcm(z1.order.clone(), z2.order.clone());
        // TODO: Use Rob's code to actually do some reductions here?
        // TODO: Reduce in-place?
        (
            z1.increase_order_to(new_order),
            z2.increase_order_to(new_order),
        )
    }
}

impl PartialEq for Cyclotomic {
    fn eq(&self, other: &Self) -> bool {
        let (new_self, new_other) = Cyclotomic::match_orders(self, other);

        // Since we assume the exponent list is sorted, we can just do a
        // comparison of the coeff and exp vectors, the number is in
        // canonical form
        new_self.coeffs == new_other.coeffs && new_self.exps == new_other.exps
    }
}

impl Eq for Cyclotomic {}

impl Add for Cyclotomic {
    type Output = Cyclotomic;

    fn add(self, other: Cyclotomic) -> Cyclotomic {
        let (new_self, new_other) = Cyclotomic::match_orders(&self, &other);

        let mut self_i = 0;
        let mut other_i = 0;
        let mut result = Cyclotomic::zero();

        // We are producing a result cyclotomic in canonical form assuming the
        // inputs are in canonical form. That is, we build the result term by
        // term, assuming that the exponent lists are sorted.

        while self_i < new_self.exps.len() || other_i < new_other.exps.len() {
            while self_i < new_self.exps.len() && new_self.exps[self_i] < new_other.exps[other_i] {
                result.coeffs.push(new_self.coeffs[self_i].clone());
                result.exps.push(new_self.exps[self_i].clone());
                self_i += 1;
            }

            while other_i < new_other.exps.len() && new_other.exps[other_i] < new_self.exps[self_i]
            {
                result.coeffs.push(new_other.coeffs[other_i].clone());
                result.exps.push(new_other.exps[other_i].clone());
                other_i += 1;
            }

            // If the exponents match up, we just add the coefficients - they are for the same term
            if new_self.exps[self_i] == new_other.exps[other_i] {
                result
                    .coeffs
                    .push(new_self.coeffs[self_i].clone() + new_other.coeffs[other_i].clone());
                result.exps.push(new_self.exps[self_i]);
                self_i += 1;
                other_i += 1
            }
        }

        result
    }
}

// \section{Improving performance}
//
// Insert our incredible mathematical, algorithmic genius here.

// \end{document}
