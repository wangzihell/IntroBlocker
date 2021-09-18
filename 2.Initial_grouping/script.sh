#!/bin/bash
set -euxo pipefail

WD=/data2/rawdata2/tetraintro/

for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B;do
  n_bin=$(wc -l < ../../bed_file2/${CHR}.5M.bed)
  
  # extract samples need to be analyzed
  for i in /data2/rawdata2/variant_density/6.cross_sample_5M_201224/raw/${CHR}*5M.dist.txt;do
    name=$(basename $i|cut -f 1,2 -d".")
    python ${WD}/bin/extra_dist.py -s order.txt -d ${i} -o ${name}.dist.txt
  done
  
  # initial grouping using hiararchical clustering algorithm, removing CNV blocks
  for (( i=1;i<=${n_bin};i+=1 ));do
    gawk -vn=$i 'NR==n{print}' /data2/rawdata2/variant_density/CNV_masker/201224/${CHR}.CNVfilter.5M.txt > ${CHR}_filter.tmp
    Rscript ${WD}/bin/do_cluster_average.R -i ${CHR}.${i}.dist.txt -f ${CHR}_filter.tmp -o ${CHR}.${i}.label.txt
  done
  
  # format output files
  rm -f ${CHR}.allsample.cluslabel
  for (( i=1;i<=${n_bin};i+=1 ));do
    datamash transpose < ${CHR}.${i}.label.txt
  done | gawk -vOFS="\t" '{for(i=1;i<=NF;i++){$i=$i+1};print}' >> ${CHR}.allsample.cluslabel
  
  ((printf "CHROM\nstart\nend\n"; cat <(cut -f 1 name.txt))| datamash transpose; ( paste ${WD}/bed_file2/${CHR}.5M.bed ${CHR}.allsample.cluslabel) ) > ${CHR}.allsample.cluslabel.fomt

done
