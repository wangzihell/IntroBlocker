# IntroBlocker: the ancestry inference algorithm
Inference and visualization of the genome ancestry mosaics leveraging the stratificaiton pattern of genomic variant density.

## Key functionalities:
1. binwise genetic distance
2. initial grouping
3. ancestry inferenc
4. bayesian smoothing
5. mosaic graph visualization

+ Run pipeline: the main functions could be evoked by the **pipeline.sh** bash script.

## Applications:
Applications of the IntroBlocker algorithm on the whole-genome resequencing dataset containing 387 interploidy wheat samples and further analysis built upon the inference could be found in the github repo [CAU-MosaicWheat](https://github.com/wangzihell/CAU-MosaicWheat).

## Demonstrations
A small demo dataset was provided in in directory Demo_data. 

## Documentation:
The documentation of IntroBlocker could be accessed at [the wiki page](https://github.com/wangzihell/IntroBlocker/wiki).  

## Prerequisites & version tested:
+ Linux system with a 64 bit CPU
+ Bcftools v1.11
+ Bedtools v2.29.2
+ Datamash v1.4
+ Plink v1.9
+ python v3.7.6
+ python packages:
  + optparse
  + numpy
+ R v3.6
+ R packages:
  + optparse
  + dendextend

## Installation & required time
clone this git repository to your local computer, and no further installation & compilation time was needed.
  
## Citation
Wang, Z., Wang, W., Xie, X., Wang, Y., Yang, Z., Peng, H., Xin, M., Yao, Y., Hu, Z., Liu, J., Su, Z., Xie, C., Li, B., Ni, Z., Sun, Q., and Guo, W. (2022). Dispersed emergence and protracted domestication of polyploid wheat uncovered by mosaic ancestral haploblock inference. Nature Communications 13: 3891.
