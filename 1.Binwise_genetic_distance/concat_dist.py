#!/usr/bin/env python

from optparse import OptionParser
import numpy as np

def concat_dist(CHR="chr1A", BINSIZE=5, THR=1000, OUTFILE="test"):
    #
    BED_PATH="/data2/rawdata2/tetraintro/bed_file2/"
    DIST_PATH="/data2/rawdata2/variant_density/5.cross_sample_201224/raw/"
    # CNV_PATH="/data2/rawdata2/variant_density/CNV_masker/"

    n_bin=sum(1 for line in open(BED_PATH+CHR+'.1M.bed'))
    n_sample=sum(1 for line in open(DIST_PATH+CHR+'.1.rawdist.dist'))
    init_mat=np.zeros((n_sample, n_sample))
    # remove CNV filter from this step, and move it to clustering step.
    # cnv_mat=np.loadtxt(CNV_PATH+CHR+".CNVfilter.txt")
    count=1
    for i in range(1, n_bin+1):  #[)
        tmp=np.loadtxt(DIST_PATH+CHR+"."+str(i)+'.rawdist.dist')
        tmp=np.nan_to_num(tmp)
        # cnv_mat_tmp=cnv_mat[i-1,]
        # tmp=cnv_mat_tmp[:,None]*cnv_mat_tmp[None,:]*tmp
        # this method dont adjust missing bin.
        init_mat=init_mat+tmp
        if i%BINSIZE==0:
            init_mat=init_mat/BINSIZE/2 # "/2: plink allele count"
            np.savetxt(CHR+"."+str(count)+".5M.dist.txt", init_mat, "%i")
            init_mat=np.zeros((n_sample, n_sample))
            count+=1
    if i%BINSIZE!=0:
        init_mat=init_mat/(i%BINSIZE)
        np.savetxt(CHR+"."+str(count)+".5M.dist.txt", init_mat, "%i")


def main():
    parser = OptionParser()
    parser.add_option("-c", dest="CHR")
    parser.add_option("-b", dest="BINSIZE", default=5)
    parser.add_option("-t", dest="THR", default=1000)
    parser.add_option("-o", dest="OUTFILE")
    (options, args) = parser.parse_args()
    #
    concat_dist(options.CHR, options.BINSIZE, options.THR, options.OUTFILE)
    
# ==========================================
if __name__ == "__main__":
    main()
