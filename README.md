# IntroBlocker: the ancestry inference algorithm
Inference and visualization of the genome ancestry mosaics leveraging the stratificaiton pattern of genomic variant density.

## Key functionalities:
1. binwise genetic distance
2. initial grouping
3. ancestry inferenc
4. bayesian smoothing
5. mosaic graph visualization

+ Run pipeline: each of the main function could be evoked by the **script.sh** bash script.

## Applications:
Applications of the IntroBlocker algorithm on the whole-genome resequencing dataset containing 387 interploidy wheat samples and further analysis built upon the inference could be found in the github repo [CAU-MosaicWheat](https://github.com/wangzihell/CAU-MosaicWheat).

## Demonstrations
A small demo dataset, corresponding instructions to run on data and expected output were provided in in directory Demo_data. 

## Documentation:
The documentation of IntroBlocker could be accessed at [readthedocs](https://introblocker.readthedocs.io/en/latest/).  
Brief descrptions of each sub-uncion could also e found in the READMD file in each folder and in the comments of scripts.  

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
  + R.utils
  + dendextend

## Installation & required time
clone this git repository to your local computer, and no further installation & compilation time was needed.
