# Cyclotomic

This is a library for exact operations in cyclotomic fields that is
well-tested, well-documented, and aims for speed in certain use
cases. Can be used from C, C++ or Rust.

## Why use this library?

There are other libraries for more general use-cases out
there, for example:

* [antic](https://github.com/wbhart/antic)
* [e-antic](https://github.com/videlec/e-antic)
* [Arb](http://arblib.org/)
* [FLINT](http://flintlib.org/)

They all do *many* things and have a wide scope. Our library only does
one thing, and does it *fast*: cyclotomic field operations. And not
just that, we are optimised for the specific use case of cyclotomic
linear algebra, with the intention of using SIMD operations wherever
it will make our code faster (backed up by extensive benchmarking).

There's a lot of structure common to all number fields, and it's
possible to get number field operations reasonably fast. But there is
much, much more that's specific to cyclotomic fields, and this allows
us to be faster for certain use cases (keep an eye on arXiv for some
details).

See [below for some benchmarks against antic](#benchmarks).

## Quick start

### Use cyclotomic in a Rust program

To install the latest release, install
[cargo](https://doc.rust-lang.org/cargo/getting-started/installation.html)
to get a functional Rust installation, including the cargo package
manager. Then add the `cyclotomic` crate to your `Cargo.toml`:

```toml
[dependencies]
cyclotomic = "1.*"
```

If you're using this before version 1 is released, then change that to
"0.*", but there'll be a lot of breaking changes.

### Use cyclotomic in a C or C++ program

TODO, but this is very important.

### More documentation

Take a look
[here](https://cyclotomicfields.github.io/cyclotomic/cyclotomic/) for
full API documentation, along with code examples.

## Quick start for devs

1. Install Rust using the instructions on their
[website](https://www.rust-lang.org/tools/install)
2. Check that Cargo is on your PATH by running `cargo --version` in the
terminal
3. Clone this repository using
`git clone https://github.com/CyclotomicFields/cyclotomic.git`
4. Download the prime tables by running the bash script
`./download_prime_tables.sh`. You may then optionally convert it to a CSV using
the J script `./convert_prime_tables_to_csv.ijs`, which will enable the program
to load the primes into memory faster
5. Check that all the tests pass using `cargo test`
6. Add your name to the contributors list in this README and push, to confirm
that you have push permissions
7. Check that your push has triggered a pipeline build on
[Travis](https://travis-ci.com/github/CyclotomicFields/cyclotomic/builds)
8. Run the benchmark by first installing nightly rust with
`rustup install nightly` and then running `cargo +nightly bench`.

## Building the documentation

Just run `cargo doc --no-deps` and have a look at
`target/doc/cyclotomic` with a web browser.

Coming soon: MathJax rendered LaTeX math in the documentation.

## Motivation

Most of the motivation for this project comes from representation
theory. In particular, given a finite group G, where m is the least
common multiple of the orders of the elements of G, let K be a
cyclotomic field containing the mth roots of unity.

Then using only coefficients from K, you can realise all
representations of G. So in some sense, the cyclotomic numbers are
"all you need" when it comes to the representation theory of finite
groups.

The goals of this project are (in order):

1. Implement field operations and equality for cyclotomic fields of
   arbitrary degree. This includes testing, since otherwise we don't
   know whether it works or not.

2. Make it readable, well documented, and not a complete mess to
   contribute to. Documentation may or may not include detailed and
   rigorous (in a mathematical sense) proofs of any
   performance-enhancing tricks used.

3. Profile and benchmark everything so the performance of this library
   is good enough that it can be used to solve non-trivial problems.

## Inspiration

The implementation we have borrowed most from is the one for the
[GAP computer algebra system](https://github.com/gap-system/gap/blob/master/src/cyclotom.c),
which implements specific logic for cyclotomic fields using various
basis tricks. This is closest in spirit to what we're doing, but isn't
standalone.

Another library is
[cyclotomic](https://hackage.haskell.org/package/cyclotomic) which is
the same thing but in Haskell. This is just a reimplementation of
GAP's algorithms, but I find Haskell code easier to read than C
sometimes, so this deserves a mention.

## Benchmarks

There are source code for some benchmarks against various libraries
under the `benchmarks/` directory. Currently, there are only
benchmarks for the antic library for number fields (among other
things).

To compile the benchmark:

```sh
$ cd benchmarks/antic
$ make
```

(assuming you have the e-antic and flint libraries installed (TODO:
maybe some guidance on how to install those))

To run, specify the order of the cyclotomic field and the number of
dot products to perform:

```sh
$ ./anticbench 100 120000
n = 100
num_tests = 120000
generating test data
starting benchmark
time elapsed (ms):
3586
```

The actual benchmark run is kind of arbitrary and will be documented
here more thoroughly. Basically some dot products.

#### Contributors

Rob Moore, Kaashif Hymabaccus

## Copyright

Copyright (C) 2020 Rob Moore, Kaashif Hymabaccus

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

The full license text can be found in the [LICENSE file](/LICENSE).
