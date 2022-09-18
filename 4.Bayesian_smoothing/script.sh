#!/bin/bash

source ../parameters.sh

while read CHR;do
  python ${script_dir}/4.Bayesian_smoothing/smoothing.py -c ${CHR} -d ../01-Binwise-genetic-distance/ -g ../03-Ancestry-inference/${CHR}.grp -m ../03-Ancestry-inference/${CHR}.ms -o ../01-Binwise-genetic-distance/plink_taxa_order.txt -v ${vote_dis_threshold} > ${CHR}.log
done < ../chr_list.txt
