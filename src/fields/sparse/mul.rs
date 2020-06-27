use crate::fields::{MultiplicativeGroup, CyclotomicFieldElement};
use crate::fields::sparse::*;

impl MultiplicativeGroup for Number {
    /// Multiplies term by term, not bothering to do anything interesting.
    fn mul(&mut self, rhs: &mut Self) -> &mut Self {
        let mut z1 = self;
        let mut z2 = rhs;
        Self::match_orders(z1, z2);

        let mut result = zero_order(z1.order.clone());

        // This order is almost certainly not optimal. But you know, whatever.
        // TODO: make it gooder
        result.order = z1.order;
        for (exp1, coeff1) in z1.coeffs.clone() {
            for (exp2, coeff2) in z2.coeffs.clone() {
                let new_exp = (exp1 + exp2) % z1.order.clone();
                let new_coeff = coeff1.clone() * coeff2.clone();

                // Special case: if the new exponent would be 0, since 1 is not
                // a basis element, we have to use the fact that:
                // $1 = -\sum_{i=1}^{p-1} \zeta_n^i$ to rewrite the new constant
                // term in our basis.
                if new_exp != 0 {
                    match result.coeffs.clone().get(&new_exp) {
                        Some(existing_coeff) => {
                            result.coeffs.insert(new_exp, new_coeff + existing_coeff)
                        }
                        None => result.coeffs.insert(new_exp, new_coeff),
                    };
                } else {
                    for i in 1..result.order.clone() {
                        match result.coeffs.clone().get(&i) {
                            Some(existing_coeff) => {
                                result.coeffs.insert(i, existing_coeff - new_coeff.clone())
                            }
                            None => result.coeffs.insert(i, -new_coeff.clone()),
                        };
                    }
                }
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
        println!("z = {:?}", z);

        // Let $L = \mathbb{Q}(\zeta_n), K = \mathbb{Q}$.
        // Then $L/K$ is a degree $\phi(n)$ extension.
        let n = z.order;

        // The Galois group $G = \text{Aut}(L/K)$ has order $\phi(n)$. The
        // elements are the automorphisms $\zeta_n \mapsto \zeta_n^i$ for all
        // $1 \leq i \leq n-1$ coprime to $n$.

        // This is the product except for the term for $t = \id_L$.
        let mut x = one_order(n);

        for i in 2..n {
            if are_coprime(i, n) {
                x.mul(&mut apply_automorphism(&z, i));
            }
        }
        println!("x = {:?}", x);

        // The full product:
        let q_cyc = z.clone().mul(&mut x.clone()).clone();
        println!("q_cyc = {:?}", q_cyc);

        // q_cyc is rational, so let's extract the rational bit (all of it).
        // We need to do some tricks to make it rational (it might be in a
        // bit of a weird form).

        // How do you tell if a cyclotomic is really rational?
        let q: Q = try_rational(&q_cyc).unwrap();
        println!("q = {:?}", q);

        let z_inv = x.scalar_mul(&q.inv());

        println!("z_inv = {:?}", z_inv);

        *self = z_inv.clone();

        self
    }
}
