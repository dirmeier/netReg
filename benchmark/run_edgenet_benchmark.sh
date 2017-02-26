#!/usr/bin/env bash

if [ "$(uname)" == "Darwin" ]; then
    exe=$(greadlink -f $0)
else
    exe=$(readlink -f $0)
fi

dir=$(dirname $exe)
scrpt="${dir}/edgenet_benchmark.R"

bsub -o ~/bewi/members/simondi/bsub/netReg.txt -e ~/bewi/members/simondi/bsub/ -W 120 -M 20000 -n 1 -R "rusage[mem=20000]" Rscript  $scrpt -n 10 -p 10
bsub -o ~/bewi/members/simondi/bsub/netReg.txt -e ~/bewi/members/simondi/bsub/ -W 120 -M 20000 -n 1 -R "rusage[mem=20000]" Rscript  $scrpt -n 100 -p 100
bsub -o ~/bewi/members/simondi/bsub/netReg.txt -e ~/bewi/members/simondi/bsub/ -W 120 -M 20000 -n 1 -R "rusage[mem=20000]" Rscript  $scrpt -n 1000 -p 100
bsub -o ~/bewi/members/simondi/bsub/netReg.txt -e ~/bewi/members/simondi/bsub/ -W 220 -M 20000 -n 1 -R "rusage[mem=20000]" Rscript  $scrpt -n 1000 -p 1000
bsub -o ~/bewi/members/simondi/bsub/netReg.txt -e ~/bewi/members/simondi/bsub/ -W 220 -M 20000 -n 1 -R "rusage[mem=20000]" Rscript  $scrpt -n 10000 -p 1000
bsub -o ~/bewi/members/simondi/bsub/netReg.txt -e ~/bewi/members/simondi/bsub/ -W 320 -M 20000 -n 1 -R "rusage[mem=20000]" Rscript  $scrpt -n 10000 -p 10000