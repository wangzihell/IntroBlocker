#!/bin/bash
set -euxo pipefail

WD=/data2/rawdata2/tetraintro/

for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B;do
  ~/anaconda2/envs/tetraintro/bin/python ${WD}/bin/find_haplo_block.py -i ../repr_dist_AABB/${CHR}.allsample.cluslabel.fomt -o ${CHR} -m sample_name.txt -M "semi-supervised" -C "priority"
done
