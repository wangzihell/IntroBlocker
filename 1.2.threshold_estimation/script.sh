#!/bin/bash

# use binwise genetic distances as the input. the file could be thinned to reduce runtime.
# the file dist.txt is a demo dataset
cat ../1.Binwise_genetic_distance/*.dist.txt |tr " " "\n" |gawk 'NR%100==1'> dist.txt

# fit and plot the distribution
Rscript plot_distribution.R -i dist.txt -o density_distribution.pdf
