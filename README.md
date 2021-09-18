# IntroBlocker: the ancestry inference algorithm
Inference and visualization of the genome ancestry mosaics leveraging the stratificaiton pattern of genomic variant density.

## key functionalities:
1. binwise genetic distance
2. initial grouping
3. ancestry inferenc
4. bayesian smoothing
5. mosaic graph visualization

+ Run pipeline: each of the main function could be evoked by the **script.sh** bash script.

## Demonstrations & Examples:
Applications of the IntroBlocker algorithm on the whole-genome resequencing dataset containing 387 interploidy wheat samples and further analysis built upon the inference could be found in the github repo [CAU-MosaicWheat](https://github.com/wangzihell/CAU-MosaicWheat).

## Documentation
Brief descrptions of each funcion could be found in the READMD file in each folder and in the comments of scripts.

## Prerequisites:
+ Linux system with a 64 bit CPU
+ Bcftools
+ Bedtools
+ Datamash
+ Plink **1.9**
+ python**3**
+ python packages:
  + optparse
  + numpy
+ R
+ R packages:
  + optparse
  + R.utils
  + dendextend

## Citation
