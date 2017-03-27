#!/usr/bin/env bash

if [ "$(uname)" == "Darwin" ]; then
    exe=$(greadlink -f $0)
else
    exe=$(readlink -f $0)
fi

dir=$(dirname $exe)
# 
scrpt="${dir}/edgenet_benchmark_time.R"
# 
echo "Doing $scrpt"
#
Rscript  $scrpt -n 1000 -p 99
Rscript  $scrpt -n 10 -p 10
Rscript  $scrpt -n 100 -p 100
Rscript  $scrpt -n 1000 -p 100
Rscript  $scrpt -n 1000 -p 1000
Rscript  $scrpt -n 10000 -p 1000


# bsub -o ~/bewi/members/simondi/bsub/netReg_time.txt -e ~/bewi/members/simondi/bsub/ -W 120 -M 20000 -n 1 -R "rusage[mem=20000]" Rscript  $scrpt -n 10 -p 10
# bsub -o ~/bewi/members/simondi/bsub/netReg_time.txt -e ~/bewi/members/simondi/bsub/ -W 120 -M 20000 -n 1 -R "rusage[mem=20000]" Rscript  $scrpt -n 100 -p 100
# bsub -o ~/bewi/members/simondi/bsub/netReg_time.txt -e ~/bewi/members/simondi/bsub/ -W 120 -M 20000 -n 1 -R "rusage[mem=20000]" Rscript  $scrpt -n 1000 -p 100
# bsub -o ~/bewi/members/simondi/bsub/netReg_time.txt -e ~/bewi/members/simondi/bsub/ -W 320 -M 20000 -n 1 -R "rusage[mem=20000]" Rscript  $scrpt -n 1000 -p 1000
# bsub -o ~/bewi/members/simondi/bsub/netReg_time.txt -e ~/bewi/members/simondi/bsub/ -W 320 -M 20000 -n 1 -R "rusage[mem=20000]" Rscript  $scrpt -n 10000 -p 1000
# bsub -o ~/bewi/members/simondi/bsub/netReg_time.txt -e ~/bewi/members/simondi/bsub/ -W 520 -M 20000 -n 1 -R "rusage[mem=20000]" Rscript  $scrpt -n 10000 -p 10000

# scrpt="${dir}/edgenet_benchmark_rss.R"
# 
# echo "Doing $scrpt" 
# 
# 
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 20:00 -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 10 -p 20 -q 10 -s 1
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 20:00 -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 10 -p 20 -q 10 -s 2
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 20:00 -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 10 -p 20 -q 10 -s 5
# 
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 24:00 -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 10 -p 30 -q 10 -s 1
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 24:00 -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 10 -p 30 -q 10 -s 2
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 24:00 -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 10 -p 30 -q 10 -s 5
# 
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 24:00 -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 100 -p 10 -q 10 -s 1
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 24:00 -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 100 -p 10 -q 10 -s 2
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 24:00 -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 100 -p 10 -q 10 -s 5
# 
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 120:00 -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 100 -p 100 -q 10 -s 1
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 120:00 -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 100 -p 100 -q 10 -s 2
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 120:00 -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 100 -p 100 -q 10 -s 5
# 
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 120:00  -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 100 -p 100 -q 100 -s 1
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 120:00  -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 100 -p 100 -q 100 -s 2
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 120:00  -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 100 -p 100 -q 100 -s 5
# 
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 120:00  -M 10000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 100 -p 1000 -q 100 -s 1
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 120:00  -M 10000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 100 -p 1000 -q 100 -s 2
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 120:00  -M 10000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 100 -p 1000 -q 100 -s 5
# 
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 120:00  -M 10000 -n 1 -R "rusage[mem=10000]" Rscript $scrpt -n 100 -p 100 -q 1000 -s 1
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 120:00  -M 10000 -n 1 -R "rusage[mem=10000]" Rscript $scrpt -n 100 -p 100 -q 1000 -s 2
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 120:00  -M 10000 -n 1 -R "rusage[mem=10000]" Rscript $scrpt -n 100 -p 100 -q 1000 -s 5
# 
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 120:00  -M 10000 -n 4 -R "rusage[mem=10000]" Rscript $scrpt -n 100 -p 1000 -q 100 -s 1
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 120:00  -M 10000 -n 4 -R "rusage[mem=10000]" Rscript $scrpt -n 100 -p 1000 -q 100 -s 2
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 120:00  -M 10000 -n 4 -R "rusage[mem=10000]" Rscript $scrpt -n 100 -p 1000 -q 100 -s 5
# 
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 120:00  -M 20000 -n 4 -R "rusage[mem=20000]" Rscript $scrpt -n 100 -p 1000 -q 1000 -s 1
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 120:00  -M 20000 -n 4 -R "rusage[mem=20000]" Rscript $scrpt -n 100 -p 1000 -q 1000 -s 2
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 120:00  -M 20000 -n 4 -R "rusage[mem=20000]" Rscript $scrpt -n 100 -p 1000 -q 1000 -s 5
# 
# 
# ##############
# 
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 20:00 -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 10 -p 20 -q 10 -s 1
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 20:00 -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 10 -p 20 -q 10 -s 2
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 20:00 -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 10 -p 20 -q 10 -s 5
# 
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 24:00 -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 10 -p 30 -q 10 -s 1
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 24:00 -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 10 -p 30 -q 10 -s 2
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 24:00 -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 10 -p 30 -q 10 -s 5
# 
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 24:00 -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 100 -p 10 -q 10 -s 1
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 24:00 -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 100 -p 10 -q 10 -s 2
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 24:00 -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 100 -p 10 -q 10 -s 5
# 
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 120:00 -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 100 -p 100 -q 10 -s 1
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 120:00 -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 100 -p 100 -q 10 -s 2
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 120:00 -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 100 -p 100 -q 10 -s 5
# 
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 120:00  -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 100 -p 100 -q 100 -s 1
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 120:00  -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 100 -p 100 -q 100 -s 2
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 120:00  -M 5000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 100 -p 100 -q 100 -s 5
# 
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 520:00  -M 10000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 100 -p 1000 -q 100 -s 1
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 520:00  -M 10000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 100 -p 1000 -q 100 -s 2
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 520:00  -M 10000 -n 4 -R "rusage[mem=5000]" Rscript $scrpt -n 100 -p 1000 -q 100 -s 5
# 
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 520:00  -M 10000 -n 1 -R "rusage[mem=10000]" Rscript $scrpt -n 100 -p 100 -q 1000 -s 1
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 520:00  -M 10000 -n 1 -R "rusage[mem=10000]" Rscript $scrpt -n 100 -p 100 -q 1000 -s 2
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 520:00  -M 10000 -n 1 -R "rusage[mem=10000]" Rscript $scrpt -n 100 -p 100 -q 1000 -s 5
# 
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 520:00  -M 10000 -n 4 -R "rusage[mem=10000]" Rscript $scrpt -n 100 -p 1000 -q 100 -s 1
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 520:00  -M 10000 -n 4 -R "rusage[mem=10000]" Rscript $scrpt -n 100 -p 1000 -q 100 -s 2
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 520:00  -M 10000 -n 4 -R "rusage[mem=10000]" Rscript $scrpt -n 100 -p 1000 -q 100 -s 5
# 
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 520:00  -M 20000 -n 4 -R "rusage[mem=20000]" Rscript $scrpt -n 100 -p 1000 -q 1000 -s 1
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 520:00  -M 20000 -n 4 -R "rusage[mem=20000]" Rscript $scrpt -n 100 -p 1000 -q 1000 -s 2
# bsub -o ~/bewi/members/simondi/bsub/ -e ~/bewi/members/simondi/bsub/ -W 520:00  -M 20000 -n 4 -R "rusage[mem=20000]" Rscript $scrpt -n 100 -p 1000 -q 1000 -s 5

