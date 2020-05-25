// This will be an implementation of GAP's cyclotomic field arithmetic algorithms. This will be the
// reference implementation for which the algorithms are assumed correct (implementation issues
// aside). All interesting mathematical facts and theorems probably come from comments in the GAP
// implementation, although some proofs and examples come from us.
extern crate num;

use std::vec::Vec;
type Q = num::rational::BigRational;

// A Cyclotomic represents an element of the union of all cyclotomic fields.
//
// Suppose $z \in \mathbb{Q}(\zeta_n)$ where $\zeta_n$ is an nth root of unity. Then it is possible
// to write $z = \sum_{i=0}^{n-1} q_i \zeta_n^i$ where $q_i \in \mathbb{Q}$. Note that the degree
// of the field extension $\mathbb{Q}(\zeta_n)/\mathbb{Q}$ is only $\phi(n)$, so this sum is not
// unique unless $n$ is prime. This is because the set $\{ \zeta_n^i : 0 \leq i \leq n-1 \}$ is not
// linearly independent over $\mathbb{Q}$. We will see later on that we can choose a particularly
// nice basis that minimises the amount of computation we have to do.
//
// If we let N = length of coefficients = length of exponents, then our encoding of $z$ is:
// $z = \sum_{i=0}^{N-1} coefficients[i] * \zeta_{order}^{exponents[i]}$.
//
// The choice of data structure here is intentionally the same as in GAP, since this is supposed
// to be the reference implementation.
#[derive(Debug)]
struct Cyclotomic {
    order: u64,
    coefficients: Vec<Q>,
    exponents: Vec<u64>,
}
