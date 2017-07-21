#!/usr/bin/env bash

OPTIND=1

netpath=""
output=""
datafold=""

show_help ()
{
    echo -e "\n$0 -n <netReg> -o <output_folder> -d <data_folder>\n"
}


if [ $# -eq 0 ];
then
    show_help
    exit 0
fi

while getopts "h?n: o: d:" opt; do
    case "$opt" in
    h)  show_help
        exit 0
        ;;
    n)  netpath=$OPTARG
        ;;
    o)  output=$OPTARG
        ;;
    d)  datafold=$OPTARG
        ;;
    *) show_help
       exit 0
       ;;
    esac
done

if [[ $netpath == "" ||  $output == "" || $datafold == "" ]];
then
    show_help
    exit 0
fi

for k in $(seq 1 1 10);
do
    X="$datafold/X_${k}.tsv"
    Y="$datafold/Y_${k}.tsv"
    GY="$datafold/GY.tsv"
    for i in $(seq 1 .5 10);
    do
        outlas="$output/out_lasso_lambda_${i}_fold_${k}.tsv"
        echo $netpath -d $X -r $Y -v $GY -l $i -o $outlas
        for j in $(seq 1 .5 20);
        do
            outedge="$output/out_netReg_lambda_${i}_xi_${j}_fold_${k}.tsv"
            echo $netpath -d $X -r $Y -v $GY -l $i -x $j -o $outedge
            exit
        done
    done
done
