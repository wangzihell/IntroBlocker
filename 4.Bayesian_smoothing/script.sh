#!/bin/bash

# settings
grp_data_folder=/data2/rawdata2/tetraintro/210109/data    # need chrXX.grp and chrXX.ms
dist_data_folder=/data2/rawdata2/tetraintro/210109/repr_dist_AABB   # need chrXX.dist.txt and name.txt
new_grp_output_folder=/data3/user3/wangwx/projs/tetraintro/210109/test
vote_dis_threshold=593 # 40*alpha = beta when dist=593

# main
WD=/data2/rawdata2/tetraintro/
TP=/data3/user3/wangwx/projs/tetraintro/210109

for CHR in chr1A chr1B chr2A chr2B chr3A chr3B chr4A chr4B chr5A chr5B chr6A chr6B chr7A chr7B;do
  python ${TP}/smoothing.py ${CHR} ${grp_data_folder} ${dist_data_folder} ${new_grp_output_folder} ${vote_dis_threshold} > ${new_grp_output_folder}/${CHR}.log
done
