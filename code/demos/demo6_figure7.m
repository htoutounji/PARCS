%% Evaluating PARCS for multivariate data (N=9) with 2 CPs and white Gaussian noise
% This code generates data for Figure 7, as
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

b = [0; 0; 0; 2; 2; 2; 0; 1; 2]'; % baselines
w = [1 2; 2 1; 2 -1; -2 0; 0 1; 0 -1; 0 0; 0 0; 0 0]'; % CP weights (for w0=1; see line 33)
[~,N] = size(w); % number of covariates
T = 100; % number of time steps
x = zeros(T,N); % time series
c = [20 60];
sd = 1; % s.d. of white Gaussian noise

% block-permutation bootstrap parameters
B = 10000; % number of bootstrap samples
k = 1    ; % block size

for r = 1:R
  rnd = randn(T,N); % white Gaussian noise
  for w0 = 0.1:0.1:1 % weight scaling parameter controlling signal-to-noise ratio
    
    % generate time series
    x(1:c(1),:)      = bsxfun(@plus,b,sd*rnd(1:c(1),:));
    x(c(1)+1:c(2),:) = bsxfun(@plus,b+w0*w(1,:),sd*rnd(1+c(1):c(2),:));
    x(c(2)+1:end,:)  = bsxfun(@plus,b+w0*(w(1,:)+w(2,:)),sd*rnd(1+c(2):end,:));
    
    % run PARCS_3, followed by block-permutation bootstrap
    model = parcs(x,3);
    model.c = c;
    model.b = b;
    model.w = w;
    model.w0 = w0;
    model.sd = sd;
    model = bpb4parcs(model,k,B,.05);
    
    % save PARCS model
    w0STR = [num2str(idivide(w0*10,int32(10))), 'p', num2str(mod(w0*10,10))];
    DataDirSTR = '../data/demo6_figure7/';
    if ~exist(DataDirSTR,'dir'), mkdir(DataDirSTR); end
    save([DataDirSTR , 'model_N',num2str(N),'_w0',w0STR,...
      '_',num2str(mod(r,R),'%03d'),'.mat'],'model');
  end
  fprintf('realization %d\n',r);
end

% (c) 2018 Hazem Toutounji, Dept. Theoretical Neuroscience, Central Institute of
% Mental Health, Medical Faculty Mannheim of Heidelberg University