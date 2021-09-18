#!/usr/bin/env Rscript

## libraries
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-i", "--infile"), dest = "infile", default = "chr1A.1.dist.txt",
              help="input file"),
  make_option(c("-o", "--outfile"), dest = "outfile", default = "chr1A.1.label.R.txt",
              help = "[opt] output file name."),
  make_option(c("-f", "--filterfile"), dest = "filterfile", default = "chr1A_filter.tmp",
              help = "[opt] output file name."),  
  make_option(c("-t", "--threshold"), dest = "threshold", default = 1000,
              help = "[opt] width of figure (inch). [default: 20]")
)

parser <- OptionParser(usage = "mapdrawer [options]",
                       option_list = option_list)

## check arguments
arguments <- parse_args(parser)
infile <- arguments$infile
outfile <- arguments$outfile
filterfile <- arguments$filterfile
threshold <- as.numeric(arguments$threshold)

DF1 <- read.table(infile)
n <- ncol(DF1)
filterfile <- as.numeric(read.table(filterfile))
DF1 <- DF1[filterfile==1, filterfile==1]
hc <- cutree(hclust(d = as.dist(DF1), method = "average"), h = threshold)
tmp <- as.numeric(names(hc))
for (i in 1:n){
  if (!(i %in% tmp)){
    hc <- append(hc, 0, after = i-1)
    names(hc)[i]=i
  }
}
hc <- hc-1
write.table(hc, outfile, quote = F, col.names = F, row.names = F)
