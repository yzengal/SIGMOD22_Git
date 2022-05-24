# Faster and Better Solution to Embed Lp Metrics by Tree Metrics
========================================================================


This repository stores the source code of our proposed framework called DC in our SIGMOD'22 paper, "Faster and Better Solution to Embed Lp Metrics by Tree Metrics".

[1] **Faster and Better Solution to Embed Lp Metrics by Tree Metrics.**
*Yuxiang Zeng, Yongxin Tong, Lei Chen.* SIGMOD 2022. 

If you find this work helpful in your research, please consider citing our paper and the bibtex is listed below:
```  
@inproceedings{DBLP:conf/sigmod/ZengTC22,
  author    = {Yuxiang Zeng and
			   Yongxin Tong and
               Lei Chen},
  title     = {Faster and Better Solution to Embed Lp Metrics by Tree Metrics},
  booktitle = {{SIGMOD}},
  year      = {2022},
}
```  

## SIGMOD22-HST.ppsx is our presentation slides.

## appendix.pdf is the additional appendix of our paper.

Usage of the source code
---------------

### Environment

gcc/g++ version: 4.8.5 with [Boost 1.74.0](https://www.boost.org/), Arya and Mount's [ANN Library 1.1.2](http://www.cs.umd.edu/~mount/ANN/), FLANN [FLANN Library 1.9.1](https://github.com/flann-lib/flann)

Python version: 2.7.5 with [Numpy 1.7.1](https://numpy.org/)

OS: Linux

### Compile the algorithm

##### Note that the directory paths of ANNLIB, BOOSTLIB and FLANNLIB in the Makefile need to be revised based on your own settings.

cd algorithm/experiment-real-and-scalability && make all
or cd algorithm/experiment-multi && make all

FRT*: the baseline FRT in our experiments. For example, FRT is for 2D space and FRT-3 is for 3D space.

Bar96*: the baseline Bar96 in our experiments. For example, Bar96 is for 2D space and Bar96-3 is for 3D space.

Bar98*: the baseline Bar98 in our experiments. For example, Bar98 is for 2D space and Bar98-3 is for 3D space.

DCnn*: our proposed DCnn in the experiments. For example, DCnn is for 2D space and DCnn-3 is for 3D space.

DCsam*: our proposed DCsam in the experiments. For example, DCsam is for 2D space and DCsam-3 is for 3D space.

DC: our proposed DC in the experiments. For example, DC is for 2D space.

HSF: the baseline HSF+FRT in our appendix (i.e., the experiment under the insertion/deletion scenario).

HSFdc: the extended algorithm HSF+DCsam in our appendix (i.e., the experiment under the insertion/deletion scenario).

FRT-5: the baseline HSF+FRT in our appendix (i.e., the experiment on the non-Lp metrics). In Makefile, "-D CHISQUARE" is for chi-square histogram distance and "-D HELLINGER" is for Hellinger distance.

DCsam-5: the extended algorithm HSF+DCsam in our appendix (i.e., the experiment on the non-Lp metrics). In Makefile, "-D CHISQUARE" is for chi-square histogram distance and "-D HELLINGER" is for Hellinger distance.

### Datasets

dataset/realData: real datasets used in the experiments. For example, checkinNYC.txt and checkinTKY are Foursquare datasets, and gaiaChengdu.txt and gaiaHaikou.txt are Didi datasets.

dataset/multiData: multi-dimensional synthetic datasets used in the experiments, eg, Uniform, Normal, Exponential and Skewed distributions.

dataset/synData: synthetic datasets for scalability tests used in the experiments, eg, Uniform, Normal, Exponential and Skewed distributions.

dataset/updateData: real datasets used in the appendix (i.e., the experiment under the insertion/deletion scenario).

dataset/nonlpData: synthetic datasets used in the appendix (i.e., the experiment on the non-Lp metrics).

genMultiData.py: the data generator of the multi-dimensional datasets

genSynData.py: the data generator of the synthetic datatests for scalability tests

genRealData.py: the data generator for varying the parameters of the real datasets


##### Notice we only provide samples here since the data file is too large. The data generators above can be used to generate all the test cases.

### Run the algorithm

**1. Single test** 

FRT ./realData/checkinNYC.txt ./realData/2/data_00.txt

HSF ./updateData/checkinNYC_2/data_0.txt

**2. Batch test**   

python run-scripts/batchRun*.py 




