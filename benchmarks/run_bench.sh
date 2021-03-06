#!/bin/bash -x

num_tests=$1
min_order=$2
max_order=$3
step=$4
d="results"
mkdir -p "$d"

for ((n=min_order; n <= max_order; n += step)); do
    echo "n=$n"
    for impl in antic sparse big_sparse; do
        echo "$n,$(../target/release/cyclobench -i $impl random -l $n -u $((n+1)) -n $num_tests -g test.g)" >> $d/$(echo $impl | tr -d '_')
    done
    echo "$n,$(gap -q < test.g)" >> $d/gap
done
