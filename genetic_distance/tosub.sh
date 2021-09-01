#!/bin/bash

for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B;do
  n=$(wc -l < ${CHR}.1M.bed)
  parallel -j procfile bash ./script.sh ${CHR} {} ::: $(eval echo {1..$n})
done
