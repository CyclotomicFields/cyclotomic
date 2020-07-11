#!/bin/bash -x

num_tests=$1
d="results"
mkdir -p "$d"

for n in {10..200..10}; do
    echo "n=$n"
    echo "$n,$(antic/anticbench $n $num_tests)" >> $d/antic
    echo "$n,$(../target/release/cyclobench -l $n -u $((n+1)) -n $num_tests -g test.g)" >> $d/cyclotomic
    echo "$n,$(gap -q < test.g)" >> $d/gap
done
