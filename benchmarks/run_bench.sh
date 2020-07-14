#!/bin/bash -x

num_tests=$1
min_order=$2
max_order=$3
step=$4
d="results"
mkdir -p "$d"

for ((n=min_order; n <= max_order; n += step)); do
    echo "n=$n"
    echo "$n,$(../target/release/cyclobench -i antic -l $n -u $((n+1)) -n $num_tests)" >> $d/antic
    echo "$n,$(../target/release/cyclobench -l $n -u $((n+1)) -n $num_tests -g test.g)" >> $d/cyclotomic
    echo "$n,$(gap -q < test.g)" >> $d/gap
done
