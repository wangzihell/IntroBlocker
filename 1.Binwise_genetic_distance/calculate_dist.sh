#!/bin/bash

source ../parameters.sh

declare -A chr_vcf
echo $chr_vcf_file
while read chr pth;do
  chr_vcf[${chr}]=${pth}
done < ${chr_vcf_file}

while read CHR;do

  num_bin=$(wc -l < ../bed/${CHR}.bed)

  for ((i=1;i<=num_bin;i++));do
    tail -n +${i} ../bed/${CHR}.bed|head -n 1 > tmp.bed
    bcftools view -v snps --min-ac=1 -M2 -m2 -R tmp.bed ${chr_vcf[${CHR}]} > tmp.vcf
    plink --vcf tmp.vcf \
        --allow-extra-chr \
        --distance square flat-missing \
        --out ${CHR}.${i}.rawdist > ${CHR}.${i}.log
  
    # convert allel count to distance
    gawk -vbin_size=${bin_size} -vOFS="\t" '{for(i=1;i<=NF;i++){$i=int($i/bin_size/2)};print}' ${CHR}.${i}.rawdist.dist > ${CHR}.${i}.dist.txt
  done

done < ../chr_list.txt

cat $(ls *.rawdist.dist.id|head -n 1)|cut -f 1 > plink_taxa_order.txt
