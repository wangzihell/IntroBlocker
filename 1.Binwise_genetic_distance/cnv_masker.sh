#!/bin/bash

source ../parameters.sh
#

declare -A taxa_bam

while read taxa chr pth;do
  taxa_bam[${taxa}_${chr}]=${pth}
done < ${taxa_bam_file}
# 

while read TAXA;do

  mkdir -p ${TAXA}

  while read CHR;do
  
    bedtools coverage -a ../bed/${CHR}.bed -b ${taxa_bam[${TAXA}_${CHR}]} -counts -sorted > ${TAXA}/${CHR}.DP
  
  done < ../chr_list.txt
  
  # normlize DP coverage by sample
  sum=$(gawk '{sum+=$4} END{print sum}' ${TAXA}/*.DP)
  line=$(wc -l ${TAXA}/*.DP | tail -n 1 | gawk '{print $1}')
  ave=$(echo "scale=4;$sum/$line"|bc)

  while read CHR;do
    gawk -vOFS="\t" -vave=$ave '{print $4/ave}' ${TAXA}/${CHR}.DP > ${TAXA}/${CHR}.norm
    sed -i '1i '$TAXA'' ${TAXA}/${CHR}.norm
  done < ../chr_list.txt

done < ../01-Binwise-genetic-distance/plink_taxa_order.txt
#

while read CHR;do
while read TAXA;do
  gawk 'NR>1{if($1<0.5 || $1>1.5){print 0}else{print 1}}' ${TAXA}/${CHR}.norm | datamash transpose
done < ../01-Binwise-genetic-distance/plink_taxa_order.txt | datamash transpose > ${CHR}.CNVfilter.txt
done < ../chr_list.txt
