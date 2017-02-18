#!/usr/bin/env bash

if [ "$(uname)" == "Darwin" ]; then
    exe=$(greadlink -f $0)
else
    exe=$(readlink -f $0)
fi

dir=$(dirname $exe)
scrpt="${dir}/edgenet_benchmark.R"

for i in 10 20
do
  Rscript $scrpt $i
done
