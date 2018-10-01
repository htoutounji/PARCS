# Detecting Multiple Change Points Using Paired Adaptive Regressors for Cumulative Sum (PARCS)

Copyright: © 2018 Hazem Toutounji.

This package is distributed under the terms of the GNU GPLv3 & Creative Commons Attribution License. Please credit the source and cite the reference above when using the code in any from of publication.

This repository contains MATLAB code for PARCS, together with demos, as described in:

[Toutounji H & Durstewitz D (2018). *Detecting multiple change points using adaptive regression splines with application to neural recordings*. Front Neuroinform. doi: 10.3389/fninf.2018.00067.](https://www.frontiersin.org/articles/10.3389/fninf.2018.00067/abstract)

```
code
└─ parcs                     %  Method Implementation
|  └─ parcs.m                %% PARCS Model Estimation
|  └─ bpb4parcs.m            %% Block-Permutation Bootstrap for PARCS
|  
└─ demos                     %  Examples
   └─ demo1_figure2.m        %% Comparing CUSUM and PARCS for white Gaussian noise
   └─ demo2_figure3.m        %% Comparing CUSUM and PARCS for moving average noise
   └─ demo3_figure4_table1.m %% Comparing maximum likelihood CUSUM and PARCS for white Gaussian noise
   └─ demo4_figure5_table2.m %% Comparing binary segmentation and PARCS for white Gaussian noise
   └─ demo5_figure6.m        %% Evaluating PARCS and block-size specification for moving average noise
   └─ demo6_figure7.m        %% Evaluating PARCS for multivariate data with white Gaussian noise
   └─ demo7_figure8.m        %% Evaluating PARCS for Poisson count data
```
