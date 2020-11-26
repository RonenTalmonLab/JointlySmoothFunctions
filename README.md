# JointlySmoothFunctions

The code in this repository creates the figures presented in this article:
https://arxiv.org/abs/2004.04386

The code was developed and tested in Matlab R2019a.
## Overview:
* Run 'MainFigure#.m' in order to generate the # figure in the paper.

## Toy problem:
* MainFigure2.m contains the implementations of Algorithm 1 and 3 (on a toy problem).
* MainFigure3.m contains the implementations of how to set the value of M.

## Sleep data:
MainFigure4.m applies Algorithm 2 to the Sleep data. It requires the precomputed kernels.

So, before running MainFigure4.m, you need to do the following:
1. Download the data from: https://physionet.org/content/capslpdb/1.0.0/
2. Download the annotation files from: https://archive.physionet.org/cgi-bin/atm/ATM
   Use the attached picture to define the correct input parameters.
3. In the file MainSleep1.m, update the 'dirPath' parameter according to where you saved the data.
4. Run MainSleep1.m, located in 'Sleep' folder. It could take several hours for it to finish running.

![Get Annotations](https://user-images.githubusercontent.com/74972592/100358313-b16f6b80-2ffe-11eb-9775-05d97408aedc.PNG)

## Non-linear dynamical systems:
* MainFigure6a.m generates the phase diagram.
  One can set the values of p1, p2, and p3 in the first cell (lines 19-21).
* MainFigure6b.m applies Algorithm 1 to the non-linear dynamical system presented in the paper.
