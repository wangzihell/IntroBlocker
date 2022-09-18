#!/usr/bin/env Rscript

## libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dendextend))

# Arguments
option_list <- list(
    make_option(c("-i", "--infile"), dest = "infile", default = NULL,
                            help = "the .grp file generated in step 3 or step 4."),
    make_option(c("-o", "--outfile"), dest = "outfile", default = NULL,
                            help = "output file name."),
    make_option(c("-b", "--binsize"), dest = "binsize", default = 5,
                            help = "the bin size in Mb [default: 5]"),
    make_option(c("-s", "--main_sample"), dest = "main_sample", default = NULL,
                            help = "[opt] main sample list to mark, like -S1,S2,S3-, without \"-\". [default: ]"), 
    make_option(c("-k", "--keepC"), dest = "keepC", default = F,
                            help = "keep color"),
    make_option(c("-L", "--localR"), dest = "localR", default = NULL,
                            help = "local region, for gene flanking."),
    make_option(c("-n", "--ncolor"), dest = "ncolor", default = 20,
                            help = "[opt] main sample list to mark, like -S1,S2,S3-, without \"-\". [default: ]"),
    make_option(c("-c", "--colorF"), dest = "colorF", default = "20_distinct_colors.txt",
                            help = "[opt] main sample list to mark, like -S1,S2,S3-, without \"-\". [default: ]"),
    make_option(c("-O", "--orderfile"), dest = "orderfile", default = NULL,
                            help = "given a sample order file."),
    make_option(c("-W", "--width"), dest = "figure.width", default = 10,
                            help = "[opt] width of figure (inch). [default: 20]"),
    make_option(c("-H", "--height"), dest = "figure.height", default = 15,
                            help = "[opt] height of figure (inch). [default: 40]")
)

parser <- OptionParser(usage = "Rscript draw_AHG.R [options]", option_list = option_list)

## check arguments
arguments <- parse_args(parser)

infile <- arguments$infile

if(infile == ""){ # default, STDIN
    infile <- file("stdin")
} else { # user specified
    if( file.access(infile) == -1){ # file not exists
        print_help(parser)
    }
}

outfile <- arguments$outfile
binsize <- as.numeric(arguments$binsize)
figure.width <- arguments$figure.width
figure.height <- arguments$figure.height
figure.format <- "pdf"
figure.resolution <- 300

## figure device
pdf(outfile, height = figure.height, width = figure.width)

#initlize graphic 
par(mar=c(0, 0, 0, 0), oma=c(1,1,1,1))

data <- read.table(file = infile, header = T, stringsAsFactors = F, check.names = FALSE, sep="\t")
left_most <- data[1,2]/binsize/1e6
right_most <- data[nrow(data),3]/binsize/1e6
data <- data[,-c(1,2,3)]

#Change ylim to change the hight of clust tree
row.height <- 0.9
upper_padding <- 5

orderfile <- arguments$orderfile
if (!is.null(orderfile)){
  sample_order <- read.table(orderfile, stringsAsFactors = F)[,1]
  tmp_index <- match(sample_order, names(data))
  # print(tmp_index)
} else {
  layout(matrix(1:2, nrow=1), widths=c(1,7))
  #clust
  clust_data <- data.frame()
  for(i in 1:ncol(data)){
      v <- c()
      for(j in 1:ncol(data)){
          sam1 <- data[,i]
          sam2 <- data[,j]
          count <- 0
          for(k in 1:length(sam1)){
              if(sam1[k] != sam2[k]){
                  count <- count+1
              }
          }
          v <- c(v,count/length(sam1))
      }
      clust_data <-rbind(clust_data, v) 
  }
  names(clust_data) <- names(data)
  hc <- hclust(d = as.dist(clust_data), method = "ward.D2")
  den <- as.dendrogram(hc)
  tmp_index <- order.dendrogram(den)
  plot(den, horiz=T, leaflab="none", yaxt = "none", dLeaf = 0, yaxs="i",ylim=c(0, ncol(data)*(row.height + 0.1)-1+upper_padding), frame.plot=F)
}

data <- data[, tmp_index]

axis.y = ncol(data)*(row.height + 0.1) + 3
nticks <- ceiling(nrow(data)*binsize/100)

bodyIndent <- 0

if ( !is.null(arguments$localR)){
  label_posi <- 5
} else {
  label_posi <- 0
}

plot(x=0, type="n", bty="n", yaxt="n", xaxt="n",
    xlab="", ylab="", 
    xlim=c(-1, nrow(data)*1.3 + bodyIndent + label_posi), ylim=c(1, ncol(data)*(row.height + 0.1)+upper_padding),
    xaxs="i", yaxs="i", frame.plot=F)

# x-axis
segments(x0 = 0 + bodyIndent, 
         y0 = axis.y, 
         x1 = nrow(data) + bodyIndent, 
         y1 = axis.y,
         col = "black",
         xpd = T
)

if ( !is.null(arguments$localR)){
  # axis ticks
  segments(x0 = c(0, nrow(data)) + bodyIndent,
           x1 = c(0, nrow(data)) + bodyIndent,
           y0 = rep(axis.y, nticks + 1),
           y1 = rep(axis.y - 0.1, nticks + 1),
           xpd = T
  )
  # axis text
  text(x = c(0, nrow(data)) + bodyIndent,
       y = rep(axis.y - 0.8, nticks + 1),
       c(left_most*binsize, right_most*binsize),
       cex = 1,
       xpd = T
  )
} else {
  # axis ticks
  segments(x0 = c(seq(from=0, to=(nticks-1) * 100, length.out = nticks), nrow(data)) + bodyIndent,
           x1 = c(seq(from=0, to=(nticks-1) * 100, length.out = nticks), nrow(data)) + bodyIndent,
           y0 = rep(axis.y, nticks + 1),
           y1 = rep(axis.y - 0.1, nticks + 1),
           xpd = T
  )
  # axis text
  text(x = c(seq(from=0, to=(nticks-1)*100, length.out = nticks)/binsize, nrow(data)) + bodyIndent,
       y = rep(axis.y - 0.8, nticks + 1),
       c(seq(from=0, to=(nticks-1)*100, length.out = nticks), nrow(data)*binsize) + left_most,
       cex = 1,
       xpd = T
  )
  
}
axis.y <- 0.2
set.seed(1)
ncolor <- as.numeric(arguments$ncolor)
keepC <- arguments$keepC

if (!is.null(arguments$colorF)){
  tmp_color_pad <- read.table(arguments$colorF, comment.char = ",", stringsAsFactors = F)[1:ncolor,1]
} else {
  tmp_color_pad <- rainbow(ncolor)[sample(ncolor)]
}
color_pad <- c("#000000", "#e5e5e5", tmp_color_pad)

data <- -abs(data)

if (arguments$keepC == F){
  groups_present <- rev(sort(unique(unlist(data))))
  if (length(groups_present) > ncolor){
    color_pad <- c(color_pad, rep("#e5e5e5", length(groups_present)-ncolor))
    # mask more than ncolor to grey
    data[data < groups_present[ncolor+1]] = 1
    names(color_pad) <- unique(c("0", "1", as.character(groups_present)))
  } else {
    data[data < groups_present[ncolor+1]] = 1
    # names(color_pad) <- unique(c("0", "1", as.character(-c(1:ncolor))))
    names(color_pad) <- unique(c("0", "1", as.character(groups_present)))
  }
} else {
  color_pad <- c(color_pad, rep("#e5e5e5", 500)) # random number , should be number of main samples
  names(color_pad) <- c("0", "1", as.character(-c(1:(500+ncolor))))
}

main_sample <- read.table(arguments$main_sample, stringsAsFactors=F)
main_sample_order <- main_sample[,2]
names(main_sample_order) <- main_sample[,1]
main_sample <- main_sample[,1]

rightPanel.extension <- 1
for(i in 1:ncol(data)){
	y.i <- i
    
	rect(xleft = 1:nrow(data) - 1 + bodyIndent, 
	     ybottom = rep(y.i*(row.height + 0.1), nrow(data)) + 0.5, 
	     xright = 1:nrow(data) + bodyIndent, 
	     ytop = rep(y.i*(row.height + 0.1), nrow(data)) + row.height-0.1 + 0.5, 
	     col = color_pad[as.character(data[,i])], 
	     border = color_pad[as.character(data[,i])])
	
    if(length(main_sample) != 0){
        for(j in 1:length(main_sample)){
            if(main_sample[j] == names(data)[i]){
                rect(xleft = mean(as.numeric(data[,i]), na.rm = T) * nrow(data)*(rightPanel.extension-1) + nrow(data) * 1.01 + 1 + bodyIndent,
                    ybottom = y.i*(row.height + 0.1) + 0.5,
                    xright = mean(as.numeric(data[,i]), na.rm = T) * nrow(data)*(rightPanel.extension-1) + nrow(data) * 1.01 + 1 + 15 + bodyIndent,
                    ytop = y.i*(row.height + 0.1) + row.height - 0.1 + 0.5,
                    col = color_pad[as.character(main_sample_order[main_sample[j]])],
                    border = NA
                    )
                break
            }
        }
    }
    text(x=mean(as.numeric(data[,i]), na.rm = T) * nrow(data)*(rightPanel.extension-1) + nrow(data) * 1.01 + 1 + bodyIndent,
		y=y.i*(row.height + 0.1)+0.1+row.height/2 + 0.5,
		names(data)[i],
		cex = 0.7, adj = c(0, 0.5)
		)
}

invisible(dev.off())