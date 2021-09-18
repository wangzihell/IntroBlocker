#!/usr/bin/env Rscript

## libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(R.utils))

# Arguments
option_list <- list(
    make_option(c("-i", "--infile"), dest = "infile", default = "",
                            help="input file"),
    make_option(c("-o", "--outfile"), dest = "outfile", default = "test.pdf",
                            help = "[opt] output file name. [default: Clustmap.SysDate.pdf]"),
    make_option(c("-m", "--mapfile"), dest = "mapfile", default = "/data2/rawdata2/sample_metadata/tetra_intro/metadata_200905.txt",
                            help="sample name transfer file. Origon name in the second column, with target name in the front.\nFile should be read-able by read.table() in R."),
    make_option(c("-b", "--binsize"), dest = "binsize", default = 5,
                            help = "[opt] extra message for debug. 0 for shut. [default: 0]"),
    make_option(c("-s", "--main_sample"), dest = "main_sample", default = "",
                            help = "[opt] main sample list to mark, like -S1,S2,S3-, without \"-\". [default: ]"), 
    make_option(c("-p", "--sample_present"), dest = "sample_present", default = NULL,
                            help = "[opt] main sample list to mark, like -S1,S2,S3-, without \"-\". [default: ]"),
    make_option(c("-k", "--keepC"), dest = "keepC", default = F,
                            help = "keep color"),
    make_option(c("-L", "--localR"), dest = "localR", default = NULL,
                            help = "local region, for gene flanking."),
    make_option(c("-n", "--ncolor"), dest = "ncolor", default = 20,
                            help = "[opt] main sample list to mark, like -S1,S2,S3-, without \"-\". [default: ]"),
    make_option(c("-c", "--colorF"), dest = "colorF", default = "/data2/rawdata2/tetraintro/bin/20_distinct_colors2.txt",
                            help = "[opt] main sample list to mark, like -S1,S2,S3-, without \"-\". [default: ]"),
    make_option(c("-g", "--groupC"), dest = "groupC", default = NULL,
                            help = "centromo cluster (.centromo.cluster)"),
    make_option(c("-G", "--geneF"), dest = "geneF", default = NULL,
                            help = "gene file"),
    make_option(c("-C", "--centro"), dest = "centro", default = NULL,
                            help = "centromere file, defining range (.centromo)"),
    make_option(c("-O", "--orderfile"), dest = "orderfile", default = NULL,
                            help = "given a sample order file."),
    make_option(c("-t", "--title"), dest = "title", default = NULL,
                            help = "[opt] title of whole file. [default: ]"),
    make_option(c("-W", "--width"), dest = "figure.width", default = 10,
                            help = "[opt] width of figure (inch). [default: 20]"),
    make_option(c("-H", "--height"), dest = "figure.height", default = 15,
                            help = "[opt] height of figure (inch). [default: 40]")
)

parser <- OptionParser(usage = "mapdrawer [options]", option_list = option_list)

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

if(outfile == ""){ # default, "FragRegView.Date"
    outfile <- paste("Clustmap", Sys.Date(), sep=".")
} else { # user specified
    outfile <- gsub(".pdf$|.png$", "", outfile, perl=T)
}

binsize <- as.numeric(arguments$binsize)
figure.width <- arguments$figure.width
figure.height <- arguments$figure.height
figure.format <- "pdf"
figure.resolution <- 300

if(! figure.format %in% c("pdf", "png")){ # format not support
    print_help(parser)
} else {
    outfile <- paste(outfile, figure.format, sep = ".")
}

## figure device
if (figure.format == "png"){
    png(outfile, height = figure.height, width = figure.width, res = figure.resolution, units = "in")
} else if (figure.format == "pdf"){
    pdf(outfile, height = figure.height, width = figure.width)
}

#initlize graphic 
if( is.null(arguments$title)){
  par(mar=c(0, 0, 0, 0), oma=c(1,1,1,1))
}else{
  par(mar=c(4, 0, 1, 0), oma=c(0,4,2,3))
}

data <- read.table(file = infile, header = T, stringsAsFactors = F, check.names = FALSE, sep="\t")
left_most <- data[1,2]/binsize/1e6
right_most <- data[nrow(data),3]/binsize/1e6
data <- data[,-c(1,2,3)]

#Change ylim to change the hight of clust tree
row.height <- 0.9
if( !is.null(arguments$geneF) ){
  upper_padding <- 7
} else {
  upper_padding <- 5
}

if (!is.null(arguments$sample_present)){
    sample_present <- read.table(arguments$sample_present, header = F)
    data <- data[,names(data) %in% sample_present$V1]
}

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
mapfile <- read.table(arguments$mapfile, header = T, stringsAsFactors = F, sep="\t")
sample_name_output <- mapfile[match(names(data), mapfile$Sample),3]

axis.y = ncol(data)*(row.height + 0.1) + 3
nticks <- ceiling(nrow(data)*binsize/100)

if( !is.null(arguments$groupC) ){
  bodyIndent <- 5
} else {
  bodyIndent <- 0
}

if ( !is.null(arguments$localR)){
  label_posi <- 5
} else {
  label_posi <- 0
}

plot(x=0, type="n", bty="n", yaxt="n", xaxt="n",
    xlab="", ylab="", 
    xlim=c(-1, nrow(data)*1.3 + bodyIndent + label_posi), ylim=c(1, ncol(data)*(row.height + 0.1)+upper_padding),
    xaxs="i", yaxs="i", main=arguments$title, frame.plot=F)

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
    names(color_pad) <- unique(c("0", "1", as.character(-c(1:ncolor))))
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
		sample_name_output[i],
		cex = 0.7, adj = c(0, 0.5)
		)
}

## group info by points
if( !is.null(arguments$groupC) ){
  groupC <- read.table(arguments$groupC, stringsAsFactors = F)
  if(length(groupC[groupC$V2=="unique",2]) >0) groupC[groupC$V2=="unique",2] <- NA
  groupC$V2 <- as.numeric(as.factor(groupC$V2))
  color <- rainbow(max(groupC$V2, na.rm = T))
  
  tmp_order <- names(data)
  groupC <- groupC[match(tmp_order, groupC$V1),]
  points(x=rep(1.5, nrow(groupC)), y=(1:ncol(data))*(row.height + 0.1)+row.height/2 + 0.5,
         col = color[groupC$V2],
         xpd = T,
         cex = 2,
         pch = 16)
  
}

## plot centromo limits
if( !is.null(arguments$centro) ){
  centro <- as.numeric(strsplit(arguments$centro, ",")[[1]])/binsize/1e6
  for (i in centro) segments(i+bodyIndent, axis.y, i+bodyIndent, ncol(data)*(row.height + 0.2)-4, lwd=2.5, lty = 2)
  segments(centro[3]+bodyIndent, axis.y, centro[3]+bodyIndent, ncol(data)*(row.height + 0.2)-4, lwd=2.5, lty = 2, col="red")
}

## plot gene location
if( !is.null(arguments$geneF) ){
  geneF <- tryCatch(read.table(arguments$geneF, stringsAsFactors = F), error=function(e) NULL)
  if (!is.null(geneF)){
    geneF <- geneF[geneF[,2]/binsize/1e6 > left_most & geneF[,3]/binsize/1e6 < right_most,]
    segments(geneF[,2]/binsize/1e6 - left_most + bodyIndent, rep(axis.y, nrow(geneF)), geneF[,2]/binsize/1e6 - left_most + bodyIndent, rep(ncol(data)*(row.height + 0.1)+6, nrow(geneF)), lty = 2)
    maptools::pointLabel(geneF[,2]/binsize/1e6 - left_most + bodyIndent, rep(ncol(data)*(row.height + 0.1)+6, nrow(geneF)), geneF[,5])
  }
}
invisible(dev.off())
