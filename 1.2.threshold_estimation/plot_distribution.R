#!/usr/bin/env Rscript

## libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(mixtools))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

option_list <- list(
  make_option(c("-i", "--infile"), dest = "infile", default = "dist.txt",
              help="binwise genetic distance (count of SNPs per 1M)"),
  make_option(c("-o", "--outfile"), dest = "outfile", default = "density_distribution.pdf",
              help = "[opt] output file name.")
)

parser <- OptionParser(usage = "mapdrawer [options]",
                       option_list = option_list)

## check arguments
arguments <- parse_args(parser)

infile <- arguments$infile
outfile <- arguments$outfile

turn_label_log <- function(x) {paste0("10e-", 6-x)}
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

# using direct hist dataset
DF1<- read.table(infile)
wait <- log(DF1$V1+1,10)
wait[wait==0]=NA
wait<-wait[!is.na(wait)]
mixmdl <- normalmixEM(wait, k = 2)

p1 <- data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), bins=50, colour = "black", fill = "grey90", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "blue", lwd = 1.5) +
  scale_x_continuous(labels = c(expression(10^-5), expression(10^-4), expression(10^-3), expression(10^-2)), limits=c(1, 4)) +
  ylab("Density") + xlab("Variant density (per bp)") +
  cowplot::theme_cowplot()
ggsave(outfile, p1, width = 5, height = 2)
