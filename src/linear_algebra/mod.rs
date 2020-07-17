use crate::fields::CyclotomicFieldElement;

/// Terrible implementation of a matrix (not contiguous)
struct Matrix<T: CyclotomicFieldElement>(Vec<Vec<T>>);

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
