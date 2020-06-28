# Cyclotomic

#### The Dream

This is a library for exact arithmetic in cyclotomic fields that is
well-tested, well-documented, and possibly fast.

Have a look at the documentation
[here](https://cyclotomicfields.github.io/cyclotomic/cyclotomic/).

## Getting Started

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

## Prior art

The best implementation of cyclotomic numbers I'm aware of is the one
in the
[GAP computer algebra system](https://github.com/gap-system/gap/blob/master/src/cyclotom.c),
which supports all of the operations we're interested in. The problem
there is that it's not standalone - it comes with a huge computer
algebra system. The other problem is that it's written in C, and
doesn't seem to be very easily parallelisable without running multiple
GAP processes.

Another library I've seen is
[cyclotomic](https://hackage.haskell.org/package/cyclotomic), which is
standalone. The problem there is that it's written in Haskell, not
focused on performance, and (in general) Haskell is hard to reason
about when it comes to performance. For me, at least.

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
