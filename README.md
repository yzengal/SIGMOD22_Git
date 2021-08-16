# Code for SIGMOD submission
========================================================================

This repository stores the source code of the proposed algorithm to solve the Embedding Lp metrics by Tree metrics (ELT) problem.

## app.pdf is the supplemental material, which provides the detailed results of FRT and DCsam in Fig. 3.

Usage of the source code
---------------

### Environment

gcc/g++ version: 4.8.5 with [Boost 1.74.0](https://www.boost.org/), Arya and Mount's [ANN Library 1.1.2](http://www.cs.umd.edu/~mount/ANN/)

Python version: 2.7.5 with [Numpy 1.7.1](https://numpy.org/)

OS: Linux

### Compile the algorithm

##### Note that the directory paths of ANNLIB and BOOSTLIB in the Makefile need to be revised based on your own settings.

cd algorithm/experiment-real-and-scalability && make all
or cd algorithm/experiment-multi && make all

FRT*: the baseline FRT in our experiments. For example, FRT is for 2D space and FRT-3 is for 3D space.

Bar96*: the baseline Bar96 in our experiments. For example, Bar96 is for 2D space and Bar96-3 is for 3D space.

Bar98*: the baseline Bar98 in our experiments. For example, Bar98 is for 2D space and Bar98-3 is for 3D space.

DCnn*: our proposed DCnn in the experiments. For example, DCnn is for 2D space and DCnn-3 is for 3D space.

DCsam*: our proposed DCsam in the experiments. For example, DCsam is for 2D space and DCsam-3 is for 3D space.

DC: our proposed DC in the experiments. For example, DC is for 2D space.

### Datasets

dataset/realData: real datasets used in the experiments. For example, checkinNYC.txt and checkinTKY are Foursquare datasets, and gaiaChengdu.txt and gaiaHaikou.txt are Didi datasets.

dataset/multiData: multi-dimensional synthetic datasets used in the experiments, eg, Uniform, Normal, Exponential and Skewed distributions.

dataset/synData: synthetic datasets for scalability tests used in the experiments, eg, Uniform, Normal, Exponential and Skewed distributions.

genMultiData.py: the data generator of the multi-dimensional datasets

genSynData.py: the data generator of the synthetic datatests for scalability tests

genRealData.py: the data generator for varying the parameters of the real datasets

##### Notice we only provide samples here since the data file is too large. The data generators above can be used to generate all the test cases.

### Run the algorithm

**1. Single test** 

FRT ./realData/checkinNYC.txt ./realData/2/data_00.txt

**2. Batch test**   

python run-scripts/batchRun*.py 




