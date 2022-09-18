#!/bin/bash

############################### the parameters which users can reset ##########################################

# the index of the reference genome with absolute path
fasta_index_file="/data2/rawdata2/demo/pseudo.fasta.fai"

# chromosomes to 
chr_list_file="/data2/rawdata2/demo/chr_list.txt"

# the working dir with absolute path
working_dir="/data2/rawdata2/demo/"

# the software with absolute path
script_dir="/data2/rawdata2/demo/IntroBlocker/"

# the bin size in Mb
bin_size=0.1

# the threshold of genomic variants density (# of variants per Mb)
thr=1000

# vote threshold
vote_dis_threshold=600

# the file contains taxa names and bam files of each chromosomes
taxa_bam_file="/data2/rawdata2/demo/taxa_bam.txt"

# the file contains vcf files of each chromosomes
chr_vcf_file="/data2/rawdata2/demo/chr_vcf.txt"

# the file contains the group of each taxa and the relative order of groups
taxa_order_file="/data2/rawdata2/demo/taxa_order.txt"

# the mode of IntroBlocker (could be "un-supervised", "semi-supervised" and "supervised")
mode="semi-supervised"

############################### end of resetting parameters ##################################################
