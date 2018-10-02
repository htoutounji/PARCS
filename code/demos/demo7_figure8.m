%% Evaluating PARCS for multivariate Poisson count data (N=9) with 2 CPs
% This code generates data for Figure 8, as
% described in Toutounji and Durstewitz (2018) Front Neuroinform.
% Please refer to 10.3389/fninf.2018.00067 for details

clc
clear
close all

rng('default');

addpath(genpath('../parcs'));

% number of realizations
R = 1000;

% time series parameters

b = [1; 1; 1; 3; 3; 3; 1; 2; 3]'; % baselines
w = [1 2; 2 1; 2 -1; -2 0; 0 1; 0 -1; 0 0; 0 0; 0 0]'; % CP weights (for w0=1; see line)
[~,N] = size(w); % number of covariates
T = 100; % number of time steps
x = zeros(T,N); % time series
c = [20 60];

% block-permutation bootstrap parameters
B = 10000; % number of bootstrap samples
k = 1    ; % block size

for r = 1:5
  % generate time series
  for n = 1:N
    x(1:c(1),n)      = poissrnd(b(n)                  ,c(1)     ,1);
    x(c(1)+1:c(2),n) = poissrnd(b(n) + w(1,n)         ,c(2)-c(1),1);
    x(c(2)+1:end,n)  = poissrnd(b(n) + w(1,n) + w(2,n),T-c(2)   ,1);
  end
  
  % run PARCS_3, followed by block-permutation bootstrap
  model = parcs(x,3);
  model.c = c;
  model.b = b;
  model.w = w;
  model = bpb4parcs(model,k,B,.05);
  
  % save PARCS model
  DataDirSTR = '../data/demo7_figure8/';
  if ~exist(DataDirSTR,'dir'), mkdir(DataDirSTR); end
  save([DataDirSTR , 'model_N',num2str(N),...
    '_',num2str(mod(r,R),'%03d'),'.mat'],'model');
  fprintf('realization %d\n',r);
end

% (c) 2018 Hazem Toutounji, Dept. Theoretical Neuroscience, Central Institute of
% Mental Health, Medical Faculty Mannheim of Heidelberg University