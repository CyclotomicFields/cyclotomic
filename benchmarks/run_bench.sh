#!/bin/bash -x

num_tests=100000

for n in {10..200..10}; do
    echo "n=$n"
    echo "$n $(antic/anticbench $n $num_tests)" >> antic_results
    echo "$n $(../target/release/cyclobench -l $n -u $((n+1)) -n $num_tests)" >> cyclotomic_results
done
