#!/usr/bin/env python

import numpy as np
from optparse import OptionParser


def main():
    usage = "Usage: %prog [-i <input file>] [-o <output file>] [-s <column number of first sample>]\n"
    #
    parser = OptionParser(usage)
    parser.add_option("-d", dest="distfile",
                  help="column id for chromosome [default: %default]", metavar="INT",
                      default = 1)
    parser.add_option("-s", dest="orderfile",
                  help="column id for chromosome [default: %default]", metavar="INT",
                      default = 1)
    parser.add_option("-o", dest="outfile",
                  help="column id for chromosome [default: %default]", metavar="INT",
                      default = 1)
    #
    (options, args) = parser.parse_args()
    #
    data_matrix = np.loadtxt(options.distfile)
    data_matrix = np.nan_to_num(data_matrix)
    sample_series = np.loadtxt(options.orderfile, "i2")
    n = len(sample_series)
    repr_dist = np.empty((n, n), "i2")

    for i in range(0, n):
        for j in range(0, n):
            # print(str(sample_series[i])+'\t'+str(sample_series[j]))
            repr_dist[i, j] = data_matrix[sample_series[i]-1, sample_series[j]-1]

    np.savetxt(options.outfile, repr_dist, "%i", '\t', '\n')

# ===========================================
if __name__ == "__main__":
    main()
