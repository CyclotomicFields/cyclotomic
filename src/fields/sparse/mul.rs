use crate::fields::{MultiplicativeGroup, CyclotomicFieldElement, Q};
use crate::fields::sparse::*;
use galois::apply_automorphism;
use basis::try_rational;

impl MultiplicativeGroup for Number {
    /// Multiplies term by term, not bothering to do anything interesting.
    fn mul(&mut self, rhs: &mut Self) -> &mut Self {
        let mut z1 = self;
        let mut z2 = rhs;
        Self::match_orders(z1, z2);

        let mut result = Number::zero_order(z1.order.clone());

        // This order is almost certainly not optimal. But you know, whatever.
        // TODO: make it gooder
        result.order = z1.order;
        for (exp1, coeff1) in &z1.coeffs {
            for (exp2, coeff2) in &z2.coeffs {
                let new_exp = (exp1 + exp2) % z1.order.clone();
                let new_coeff = coeff1 * coeff2;
                add_single(&mut result.coeffs, new_exp, new_coeff.clone());
            }
        }

        *z1 = result;
        z1
    }

    /// Gives the inverse of $z$ using the product of Galois conjugates.
    ///
    /// I don't think there's a "trivial" or "stupid" way of doing this.
    /// The product of the Galois conjugates is rational, we can normalise
    /// to get the multiplicative inverse.
    fn mul_invert(&mut self) -> &mut Self {
        let z = self.clone();
        // Let $L = \mathbb{Q}(\zeta_n), K = \mathbb{Q}$.
        // Then $L/K$ is a degree $\phi(n)$ extension.
        let n = z.order;

        // The Galois group $G = \text{Aut}(L/K)$ has order $\phi(n)$. The
        // elements are the automorphisms $\zeta_n \mapsto \zeta_n^i$ for all
        // $1 \leq i \leq n-1$ coprime to $n$.

        // This is the product except for the term for $t = \id_L$.
        let mut x = Number::one_order(n);

        for i in 2..n {
            if are_coprime(i, n) {
                x.mul(&mut apply_automorphism(&z, i));
            }
        }

        // The full product:
        let q_cyc = z.clone().mul(&mut x.clone()).clone();

        // q_cyc is rational, so let's extract the rational bit (all of it).
        // We need to do some tricks to make it rational (it might be in a
        // bit of a weird form).

        // How do you tell if a cyclotomic is really rational?
        let q: Q = try_rational(&q_cyc).unwrap();

        let z_inv = x.scalar_mul(&q.inv());

        *self = z_inv.clone();

        self
    }
}
