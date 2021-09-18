#!/bin/bash
set -euxo

WD="/data2/rawdata2/mergeFile/mergeTetra/all"
BINSIZE="5M"

# make bed files for genomic windows
bedtools makewindows -g /data/genome/wheat/CS_IWGSC/v1/161010_Chinese_Spring_v1.0_pseudomolecules.fasta.fai -w 1000000 > CS_win1M.bed
gawk '{print > $1".1M.bed"}' CS_win1M.bed

for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B;do

  plink --bcf <(bcftools view -v snps --min-ac=1 -M2 -m2 -R <(tail -n +${N} ${CHR}.1M.bed|head -n 1) ${CHR}.ann.bcf.gz -Ob) \
      --allow-extra-chr \
      --distance square flat-missing \
      --out ${CHR}.${N}.rawdist

  # convert allel count to distance
  python ./concat_dist.py -c ${CHR} &

  # detect CNV blocks
  bash ./CNV_masker.sh 
done
