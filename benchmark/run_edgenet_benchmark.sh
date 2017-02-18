#!/usr/bin/env bash

if [ "$(uname)" == "Darwin" ]; then
    exe=$(greadlink -f $0)
else
    exe=$(readlink -f $0)
fi

dir=$(dirname $exe)
scrpt="${dir}/edgenet_benchmark.R"

for i in 10 20 50 100 1000 10000
do
  bsub -o ~/bewi/members/simondi/bsub/netReg.txt -e ~/bewi/members/simondi/bsub/ -W 120 -M 20000 -n 1 -R "rusage[mem=20000]" Rscript  $scrpt $i
done
