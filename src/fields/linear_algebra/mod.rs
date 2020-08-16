use crate::fields::{CyclotomicFieldElement, FieldElement};
use std::marker::PhantomData;

/// Terrible implementation of a matrix (not contiguous)
pub struct Matrix<T, E> {
    pub value: Vec<Vec<T>>,
    pub exp: PhantomData<E>,
}

// TODO: more tests! or indeed any tests at all

impl<T, E> Matrix<T, E>
where
    T: CyclotomicFieldElement<E>,
    E: From<i64>,
{
    pub fn zero_matrix(N: usize, M: usize) -> Self {
        Matrix {
            value: vec![vec![T::zero_order(E::from(1)); M]; N],
            exp: PhantomData
        }
    }

    // should this be called one_matrix for consistency? really makes you think
    pub fn identity_matrix(N: usize) -> Self {
        let mut result = Self::zero_matrix(N, N);
        for i in 0..N {
            result.value[i][i] = T::zero_order(E::from(1));
        }
        result
    }
    fn assert_compatible(mA: &Self, mB: &Self) {
        let A = &mA.value;
        let B = &mB.value;

        assert!(A.len() > 0);
        assert!(A.into_iter().all(|v| { v.len() > 0 }));
        assert!(B.len() > 0);
        assert!(B.into_iter().all(|v| { v.len() > 0 }));

        for row in A {
            assert_eq!(row.len(), B.len());
        }
    }

    pub fn mul(mA: &mut Self, mB: &mut Self) -> Self {
        Self::assert_compatible(mA, mB);
        let A = &mut mA.value;
        let B = &mut mB.value;

        let N = A.len();
        let M = B.len();
        let L = B[0].len();

        let mut result = Self::zero_matrix(N, L);

        for i in 0..N {
            for j in 0..L {
                let mut sum = T::zero_order(E::from(1));
                for k in 0..M {
                    sum.add(&mut A[i][k].clone().mul(&mut B[k][j]));
                }
                result.value[i][j] = sum;
            }
        }

        result
    }

    pub fn add(mA: &mut Self, mB: &mut Self) -> Self {
        // just assume they're compatible
        // TODO: do something better

        let A = &mut mA.value;
        let B = &mut mB.value;

        let N = A.len();
        let M = A[0].len();

        let mut result = A.clone();

        for i in 0..N {
            for j in 0..M {
                result[i][j].add(&mut B[i][j]);
            }
        }

        Matrix {
            value: result,
            exp: PhantomData
        }
    }
}

struct Vector<T, E> {
    value: Vec<T>,
    exp: PhantomData<E>,
}

impl<T, E> Vector<T, E>
where
    T: CyclotomicFieldElement<E>,
    E: From<i64>,
{
    fn dot_product(&self, other: &mut Self) -> T {
        let v1 = &self.value;
        let v2 = &mut other.value;
        let mut result = T::zero_order(E::from(1));

        assert_eq!(v1.len(), v2.len());

        for i in 0..v1.len() {
            result.add(&mut v1[i].clone().add(&mut v2[i]));
        }

        result
    }

    // NOTE: the L2 norm might not be cyclotomic, so we can't
    // have a norm function without general number fields.
}
