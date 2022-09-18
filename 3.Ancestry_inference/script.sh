#!/bin/bash

source ../parameters.sh

tail -n +2 ${taxa_order_file} > taxa_order.txt
while read CHR;do
  python ${script_dir}/3.Ancestry_inference/find_haplo_block.py -i ../02-Initial-grouping/${CHR}.allsample.cluslabel.fomt -o ${CHR} -m taxa_order.txt -M ${mode}
done < ../chr_list.txt
