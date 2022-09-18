#!/bin/bash

source parameters.sh

mkdir -p bed
cd bed

bin_size_bp=$(echo "${bin_size}*1000000"|bc)
# make bed files for genomic windows
bedtools makewindows -g ${fasta_index_file} -w ${bin_size_bp} > ref.bed

# split into chromosomes
gawk '{print > $1".bed"}' ref.bed
