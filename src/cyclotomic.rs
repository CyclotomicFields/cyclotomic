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
// If we let $N = |\text{coefficients}| = |\text{exponents}|$, then
// our encoding of $z$ is:
//
// $$z = \sum_{i=0}^{N-1} \text{coefficients}[i]
// \zeta_{\text{order}}^{\text{exponents[i]}}$$
//
// Although we are implementing this library in Rust, not C, we are
// aiming to match GAP's implementation as closely as possible, for
// comparison purposes with later encodings and algorithms.
extern crate num;

use std::vec::Vec;
type Q = num::rational::BigRational;

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
