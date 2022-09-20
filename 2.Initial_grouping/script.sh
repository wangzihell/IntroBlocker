#!/bin/bash

source ../parameters.sh

while read CHR;do

  num_bin=$(wc -l < ../bed/${CHR}.bed)
  
  # initial grouping using hiararchical clustering algorithm, removing CNV blocks
  for (( i=1;i<=${num_bin};i+=1 ));do
    gawk -vn=$i 'NR==n{print}' ../cnv_masker/${CHR}.CNVfilter.txt > ${CHR}_filter.tmp
    Rscript ${script_dir}/2.Initial_grouping/do_cluster_average.R -i ../01-Binwise-genetic-distance/${CHR}.${i}.dist.txt -t ${cluster_thr} -f ${CHR}_filter.tmp -o ${CHR}.${i}.label.txt
  done
  
  # format output files
  rm -f ${CHR}.allsample.cluslabel
  for (( i=1;i<=${num_bin};i+=1 ));do
    datamash transpose < ${CHR}.${i}.label.txt
  done | gawk -vOFS="\t" '{for(i=1;i<=NF;i++){$i=$i+1};print}' >> ${CHR}.allsample.cluslabel
  
  ((printf "CHROM\nstart\nend\n"; cat <(cut -f 1 ../01-Binwise-genetic-distance/plink_taxa_order.txt))| datamash transpose; ( paste ../bed/${CHR}.bed ${CHR}.allsample.cluslabel) ) > ${CHR}.allsample.cluslabel.fomt

done < ../chr_list.txt
