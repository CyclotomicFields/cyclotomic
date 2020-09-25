// The absolute bare bones needed to do a character inner product

use crate::fields::rational::Rational;
use crate::fields::sparse::basis::{convert_to_base, try_reduce};
use crate::fields::{
    sparse, AdditiveGroupElement, CyclotomicFieldElement, FieldElement, MultiplicativeGroupElement,
    Z,
};
use core::arch::x86_64::*;
use faster::*;
use std::alloc::{alloc, alloc_zeroed, dealloc};
use std::alloc::Layout;
use std::{mem, ptr};
use std::slice;

#[derive(Clone)]
pub struct Number {
    /// Dense vector of coefficients, no tricks
    coeffs: Vec<f32>,
}

impl Number {
    pub fn nearest_int(&mut self) -> Z {
        // This probably works at least some of the time
        let n = self.coeffs.len();
        let mut int_coeffs = sparse::ExpCoeffMap::<i64, rug::Rational>::default();

        for i in 0..n {
            let int_rat = rug::Rational::from((self.coeffs[i].round() as i64, 1));

            int_coeffs.insert(i as i64, int_rat);
        }

        let mut z = sparse::Number::<i64, rug::Rational>::new(&(n as i64), &int_coeffs);
        z = convert_to_base(&mut z);
        try_reduce(&mut z);

        for (_, coeff) in z.coeffs {
            if !coeff.is_zero() {
                return coeff.numer().clone();
            }
        }

        return Z::from(0);
    }

    pub fn increase_order_to(z: &mut Self, new_order: i64) {
        let mut new_coeffs = vec![0_f32; new_order as usize];
        let n = z.coeffs.len();
        for exp in 0..n {
            new_coeffs[new_order as usize * exp / n] = z.coeffs[exp];
        }
        z.coeffs = new_coeffs;
    }

    pub fn match_orders(z1: &mut Number, z2: &mut Number) {
        if z1.coeffs.len() == z2.coeffs.len() {
            return;
        }
        let new_order = num::integer::lcm(z1.coeffs.len(), z2.coeffs.len());
        Number::increase_order_to(z1, new_order as i64);
        Number::increase_order_to(z2, new_order as i64);
    }
}

impl AdditiveGroupElement for Number {
    fn add(&mut self, z: &mut Self) -> &mut Self {
        Number::match_orders(self, z);
        let n = self.coeffs.len();

        self.coeffs = (
            self.coeffs.as_slice().simd_iter(f32s(0_f32)),
            z.coeffs.as_slice().simd_iter(f32s(0_f32)),
        )
            .zip()
            .simd_map(|(a, b)| a + b)
            .scalar_collect();

        self
    }

    fn add_invert(&mut self) -> &mut Self {
        panic!("unimplemented")
    }
}

unsafe fn read_leftover(arr: &[f32]) -> __m256 {
    // precondition, 0 < arr.len() < 8
    match arr.len() {
        1 => _mm256_set_ps(arr[0], 0_f32, 0_f32, 0_f32, 0_f32, 0_f32, 0_f32, 0_f32),
        2 => _mm256_set_ps(arr[0], arr[1], 0_f32, 0_f32, 0_f32, 0_f32, 0_f32, 0_f32),
        3 => _mm256_set_ps(arr[0], arr[1], arr[2], 0_f32, 0_f32, 0_f32, 0_f32, 0_f32),
        4 => _mm256_set_ps(arr[0], arr[1], arr[2], arr[3], 0_f32, 0_f32, 0_f32, 0_f32),
        5 => _mm256_set_ps(arr[0], arr[1], arr[2], arr[3], arr[4], 0_f32, 0_f32, 0_f32),
        6 => _mm256_set_ps(arr[0], arr[1], arr[2], arr[3], arr[4], arr[5], 0_f32, 0_f32),
        7 => _mm256_set_ps(
            arr[0], arr[1], arr[2], arr[3], arr[4], arr[5], arr[6], 0_f32,
        ),
        _ => panic!("leftover doesnt't have size in [1..7]"),
    }
}

impl MultiplicativeGroupElement for Number {
    fn mul(&mut self, z: &mut Self) -> &mut Self {
        Number::match_orders(self, z);
        let n = self.coeffs.len();

        // This is so we can read the shifted z directly without doing
        // i % n everywhere.
        let z_twice = [z.coeffs.as_slice(), z.coeffs.as_slice()].concat();
        let num_full_chunks = n / 8;

        // yes I know assuming AVX2 is bad, do not run this code on a Pentium 3
        unsafe {
            // TODO: make everything aligned to 32 bytes including the original coeffs

            // needed because we want to operate on chunks of 8 for the whole
            // vector
            let rounded_n = ((n + 7) / 8) * 8;

            let leftover = n % 8;

            let mem_layout = Layout::from_size_align(rounded_n * mem::size_of::<f32>(), 32).unwrap();

            // Will accumulate the final result
            let result = alloc_zeroed(mem_layout.clone()) as *mut f32;

            // After each shift, shifted_result[i] contains the coefficient of
            // x^(2i + shift mod n).
            let shifted_result = alloc(mem_layout.clone()) as *mut f32;

            // We then unjumble so unshifted_result[i] contains the coefficient of
            // x^i.
            let mut unshifted_result = alloc(mem_layout.clone()) as *mut f32;

            eprintln!("loop");
            for shift in 0..n {
                let z_shifted = &z_twice[shift..shift + n];

                // shifted_result
                for i in 0..num_full_chunks {
                    eprintln!("---");
                    let self_chunk = _mm256_loadu_ps(self.coeffs[i * 8..(i + 1) * 8].as_ptr());
                    eprintln!("self_chunk: {:?}", self_chunk);
                    let z_chunk = _mm256_loadu_ps(z_shifted[i * 8..(i + 1) * 8].as_ptr());
                    eprintln!("z_chunk   : {:?}", z_chunk);
                    let prod = _mm256_mul_ps(self_chunk, z_chunk);
                    eprintln!("prod      : {:?}", prod);
                    _mm256_storeu_ps(shifted_result.add(i * 8), prod);
                }
                if leftover != 0 {
                    let last_self_chunk = read_leftover(&self.coeffs[n - leftover..n]);
                    let last_z_chunk = read_leftover(&z_shifted[n - leftover..n]);
                    let prod = _mm256_mul_ps(last_self_chunk, last_z_chunk);
                    _mm256_storeu_ps(shifted_result.add(num_full_chunks * 8), prod);
                }

                // There is a 2* here, surely there is some cool hack I can do to
                // make this faster.
                // TODO: can we really do this without initialising the memory? is this
                // a bijection?

                ptr::write_bytes::<f32>(unshifted_result, 0, rounded_n);
                for i in 0..n {
                    *unshifted_result.add((2 * i + shift) % n) = *shifted_result.add(i);
                }
                eprintln!("shifted  : {:?}", slice::from_raw_parts::<f32>(shifted_result, rounded_n));
                eprintln!("unshifted: {:?}", slice::from_raw_parts::<f32>(unshifted_result, rounded_n));

                for i in 0..num_full_chunks {
                    let result_chunk = _mm256_loadu_ps(result.add(i * 8));
                    let unshifted_result_chunk = _mm256_loadu_ps(unshifted_result.add(i * 8));
                    let sum = _mm256_add_ps(result_chunk, unshifted_result_chunk);
                    _mm256_storeu_ps(result.add(i * 8), sum);
                }
                if leftover != 0 {
                    let last_result_chunk = read_leftover(slice::from_raw_parts::<f32>(result.add(n-leftover), leftover));
                    let last_unshifted_chunk = read_leftover(slice::from_raw_parts::<f32>(unshifted_result.add(n-leftover), leftover));
                    let sum = _mm256_mul_ps(last_result_chunk, last_unshifted_chunk);
                    _mm256_storeu_ps(unshifted_result.add(num_full_chunks * 8), sum);
                }
            }

            for i in 0..n {
                self.coeffs[i] = *result.add(i);
            }
        }

        self
    }

    fn mul_invert(&mut self) -> &mut Self {
        panic!("unimplemented")
    }
}

impl FieldElement for Number {
    fn eq(&mut self, other: &mut Self) -> bool {
        panic!("unimplemented")
    }
}

impl CyclotomicFieldElement<i64, rug::Rational> for Number {
    fn e(n: &i64, k: &i64) -> Self {
        let mut coeffs = vec![0_f32; *n as usize];
        coeffs[*k as usize] = 1_f32;
        Number { coeffs: coeffs }
    }

    fn scalar_mul(&mut self, scalar: &rug::Rational) -> &mut Self {
        let q = scalar.to_f32();
        for i in 0..self.coeffs.len() {
            self.coeffs[i] *= q;
        }
        self
    }

    fn zero_order(n: &i64) -> Number {
        Number {
            coeffs: vec![0_f32; *n as usize],
        }
    }

    fn one_order(n: &i64) -> Number {
        let mut coeffs = vec![0_f32; *n as usize];
        coeffs[0] = 1_f32;
        Number { coeffs: coeffs }
    }

    fn complex_conjugate(&self) -> Self {
        let n = self.coeffs.len();
        let mut new_coeffs = vec![0_f32; n];

        new_coeffs[0] = self.coeffs[0];
        for exp in 1..n {
            new_coeffs[(n - exp) as usize] = self.coeffs[exp as usize];
        }

        Number { coeffs: new_coeffs }
    }
}
