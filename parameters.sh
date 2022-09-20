#!/bin/bash

############################### the parameters which users can reset ##########################################

# the absolute path of index file (.fai file) of the reference genome
fasta_index_file="/data2/rawdata2/demo/pseudo.fasta.fai"

# chromosomes to analyze
chr_list_file="/data2/rawdata2/demo/chr_list.txt"

# the absolute path of the working directory
working_dir="/data2/rawdata2/demo/"

# the absolute path of the IntroBlocker software directory
script_dir="/data2/rawdata2/demo/IntroBlocker/"

# the AHG bin size (in Mb)
bin_size=0.1

# the genomic variants density threshold used in the initial grouping (# of variants per Mb)
cluster_thr=1000

# the genomic variants density threshold used in the Bayesian smoothing (# of variants per Mb)
smooth_thr=600

# the taxa_bam_file file contains taxa names and absolte paths of corresponding bam files of each chromosome
taxa_bam_file="/data2/rawdata2/demo/taxa_bam.txt"

# the chr_vcf_file file contains absolte paths of vcf files of each chromosome
chr_vcf_file="/data2/rawdata2/demo/chr_vcf.txt"

# the taxa_order_file file specifies the group of each taxa
taxa_order_file="/data2/rawdata2/demo/taxa_order.txt"

# the mode of IntroBlocker (could be "un-supervised", "semi-supervised" and "supervised")
mode="semi-supervised"

############################### end of resetting parameters ##################################################
