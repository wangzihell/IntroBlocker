#!/bin/bash

CHR=$1
N=$2

plink --bcf <(bcftools view -v snps --min-ac=1 -M2 -m2 -R <(tail -n +${N} ${CHR}.1M.bed|head -n 1) ${CHR}.ann.bcf.gz -Ob) \
      --allow-extra-chr \
      --distance square flat-missing \
      --out ${CHR}.${N}.rawdist
