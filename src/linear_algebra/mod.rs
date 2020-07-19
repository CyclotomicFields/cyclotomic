use crate::fields::CyclotomicFieldElement;

/// Terrible implementation of a matrix (not contiguous)
struct Matrix<T: CyclotomicFieldElement>(Vec<Vec<T>>);

// TODO: more tests!

impl<T> Matrix<T>
where
    T: CyclotomicFieldElement,
{
    fn assert_compatible(mA: &Self, mB: &Self) {
        let A = &mA.0;
        let B = &mB.0;

        assert!(A.len() > 0);
        assert!(A.into_iter().all(|v| { v.len() > 0 }));
        assert!(B.len() > 0);
        assert!(B.into_iter().all(|v| { v.len() > 0 }));

        for row in A {
            assert_eq!(row.len(), B.len());
        }
    }

    fn mul(mA: &mut Self, mB: &mut Self) -> Self {
        Self::assert_compatible(mA, mB);
        let A = &mut mA.0;
        let B = &mut mB.0;

        let N = A.len();
        let M = B.len();
        let L = B[0].len();

        let mut result = vec![vec![T::zero_order(1); M]; N];

        for i in 0..N {
            for j in 0..L {
                let mut sum = T::zero_order(1);
                for k in 0..M {
                    sum.add(&mut A[i][k].clone().mul(&mut B[k][j]));
                }
                result[i][j] = sum;
            }
        }

        Matrix(result)
    }
}

struct Vector<T: CyclotomicFieldElement>(Vec<T>);

impl<T> Vector<T>
where
    T: CyclotomicFieldElement,
{
    fn dot_product(&self, other: &mut Self) -> T {
        let v1 = &self.0;
        let mut v2 = &mut other.0;
        let mut result = T::zero_order(1);

        assert_eq!(v1.len(), v2.len());

        for i in 0..v1.len() {
            result.add(&mut v1[i].clone().add(&mut v2[i]));
        }

        result
    }

    // NOTE: the L2 norm might not be cyclotomic, so we can't
    // have a norm function without general number fields.
}
