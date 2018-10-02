# Detecting Multiple Change Points Using Paired Adaptive Regressors for Cumulative Sum (PARCS)

This package is distributed under the terms of the GNU GPLv3 & Creative Commons Attribution License. Please credit the source and cite the reference [below](https://www.frontiersin.org/articles/10.3389/fninf.2018.00067/abstract) when using the code in any from of publication.

This repository contains MATLAB code for PARCS, together with demos, as described in:

[Toutounji H and Durstewitz D (2018) *Detecting Multiple Change Points Using Adaptive Regression Splines With Application to Neural Recordings*. Front. Neuroinform. 12:67. doi: 10.3389/fninf.2018.00067](https://www.frontiersin.org/articles/10.3389/fninf.2018.00067/abstract)

Copyright © 2018 Toutounji and Durstewitz

```
code
├─ demos                      %  examples
│  ├─ demo1_figure2.m         %   comparing CUSUM and PARCS for white Gaussian noise
│  ├─ demo2_figure3.m         %   comparing CUSUM and PARCS for moving average noise
│  ├─ demo3_figure4_table1.m  %   comparing maximum likelihood CUSUM and PARCS for white Gaussian noise
│  ├─ demo4_figure5_table2.m  %   comparing binary segmentation and PARCS for white Gaussian noise
│  ├─ demo5_figure6.m         %   evaluating PARCS and block-size specification for moving average noise
│  ├─ demo6_figure7.m         %   evaluating PARCS for multivariate data with white Gaussian noise
│  └─ demo7_figure8.m         %   evaluating PARCS for Poisson count data
└─ parcs                      %  method implementation
   ├─ parcs.m                 %   PARCS model estimation
   └─ bpb4parcs.m             %   block-permutation bootstrap for PARCS
```
