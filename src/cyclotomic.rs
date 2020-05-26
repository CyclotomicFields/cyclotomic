// \documentclass{article}
// \usepackage{amsfonts}
// \usepackage{amsmath}
// \usepackage{parskip}
// \usepackage{geometry}[margin=1in]
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
// Currently, there is no standalone, high perfomance, open source
// implementation of cyclotomic field arithmetic. Existing
// implementations are usually part of a computer algebra system and
// thus not standalone. Examples of this are GAP's and Magma's
// implementations. This makes it difficult to benchmark and compare
// various encodings of cyclotomic fields, as well as algorithms for
// the field operations. Our objective is to produce a well-tested,
// well-documented, performant library for arithmetic with the
// cyclotomic numbers.
//
// \end{abstract}
//
// \section{Introduction}
//
// We used this library for $\mathbb{Q}$ due to performance or
// something, idk.
//
extern crate num;

use std::vec::Vec;
type Q = num::rational::BigRational;

// \section{Implementation of GAP's algorithms}
//
// Suppose $z \in \mathbb{Q}(\zeta_n)$ where $\zeta_n$ is an nth root
// of unity.  Then it is possible to write $z = \sum_{i=0}^{n-1} q_i
// \zeta_n^i$ where $q_i \in \mathbb{Q}$. Note that the degree of the
// field extension $\mathbb{Q}(\zeta_n)/\mathbb{Q}$ is only $\phi(n)$,
// so this sum is not unique unless $n$ is prime. This is because the
// set $\{ \zeta_n^i : 0 \leq i \leq n-1 \}$ is not linearly
// independent over $\mathbb{Q}$. We will see later on that we can
// choose a particularly nice basis that minimises the amount of
// computation we have to do.
//
// If we let $N = \# \text{coefficients} = \# \text{exponents}$, then
// our encoding of $z$ is:
//
// $$z = \sum_{i=0}^{N-1} \text{coefficients}[i]
// \zeta_{\text{order}}^{\text{exponents[i]}}$$
//
// Although we are implementing this algorithm in Rust, not C, we are
// aiming to match GAP's algorithm as closely as possible - this is
// the reference implementation.
#[derive(Debug)]
struct Cyclotomic {
    order: u64,
    coefficients: Vec<Q>,
    exponents: Vec<u64>,
}

// \section{Improving performance}
//
// Insert our incredible mathematical, algorithmic genius here.

// \end{document}
