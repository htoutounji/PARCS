%% Comparing CUSUM_ML and PARCS for AMOC and white Gaussian noise
% This code generates data for Figure 4 and Table 1, as
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
b = 0; % baseline
w = 1; % CP weight
sig = 1; % s.d. of white Gaussian noise

% block-permutation bootstrap parameters
B = 10000; % number of bootstrap samples
k = 1    ; % block size

% for block-permutation bootstrap of CUSUM
SbsCUSUM = zeros(B,1);

for r = 1:R
  for T = [26 , 50 , 100] % number of time steps
    
    p = sqrt(T./((1:T-1).*(T-(1:T-1))))';
    N   = (T-1)/k; % number of blocks
    prmMat = reshape(1:T-1,k,N);
    
    for cp = round((0.2:.1:.8)*T)
      
      % generate time series
      x = zeros(T,1);
      x(1:cp) = b + sig*randn(cp,1);
      x(cp+1:end) = b + w + sig*randn(T-cp,1);
      
      % run PARCS_1, followed by block-permutation bootstrap
      model = parcs(x,1);
      model.c = cp;
      model.b = b;
      model.w = w;
      model.sd = sig;
      model = bpb4parcs(model,k,B,.05);
      
      % CUSUM_ML-transformed time series
      x = model.x(1:end-1);
      y = p.*cumsum(x-mean(x));
      
      % bootstrap CUSUM_ML
      [Scs,chCS] = max(abs(y));
      bh = mean(x(1:chCS));
      wh = mean(x(chCS+1:end)) - bh;
      xST = [x(1:chCS); x(chCS+1:end) - wh];
      for bs = 1:B
        bsInd = prmMat(:,randperm(N));
        xSTbsCUSUM = xST(bsInd);
        xSTbsCUSUM = xSTbsCUSUM(1:T-1);
        ySTbsCUSUM = p.*cumsum(xSTbsCUSUM-mean(xSTbsCUSUM));
        SbsCUSUM(bs) = max(abs(ySTbsCUSUM));
      end
      
      % test for CUSUM_ML's CP significance
      model.chCS   = chCS;
      model.chCSBS = [];
      if Scs > quantile(SbsCUSUM,1-.05)
        model.chCSBS = chCS;
      end
      
      % save model (includes both CUSUM_ML and PARCS results)
      DataDirSTR = '../data/demo3_figure4_table1/';
      if ~exist(DataDirSTR,'dir'), mkdir(DataDirSTR); end
      save([DataDirSTR , 'model_T',num2str(T,'%03d'),'_c',num2str(cp,'%02d'),...
        '_',num2str(mod(r,R),'%03d'),'.mat'],'model');
      
    end
  end
  fprintf('realization %d\n',r);
end

% (c) 2018 Hazem Toutounji, Dept. Theoretical Neuroscience, Central Institute of
% Mental Health, Medical Faculty Mannheim of Heidelberg University