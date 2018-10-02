%% Comparing CUSUM and PARCS for AMOC and white Gaussian noise
% This code generates data for Figure 2, as
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
T = 100; % number of time steps
x = zeros(T,1); % time series
Cs = 20:10:80; C = numel(Cs); % ground truth CPs
Ss = .4:.1:1 ; S = numel(Ss); % s.d. of white Gaussian noise

% block-permutation bootstrap parameters
B = 10000; % number of bootstrap samples
k = 1    ; % block size
N = T/k  ; % number of blocks

% for block-permutation bootstrap of CUSUM
prmMat = reshape(1:T,k,N);
SbsCUSUM = zeros(B,1);

for r = 1:R
  for c = 1:C
    for s = 1:S
      sd = Ss(s);
      cp  = Cs(c);
      
      % generate time series
      x(1:cp) = b + sd*randn(cp,1);
      x(cp+1:end) = b + w + sd*randn(T-cp,1);
      
      % run PARCS_1, followed by block-permutation bootstrap
      model = parcs(x,1);
      model.c = cp;
      model.b = b;
      model.w = w;
      model.sd = sd;
      model = bpb4parcs(model,k,B,.05);
      
      % CUSUM-transformed time series
      y = model.y;
      
      % bootstrap CUSUM
      [Scs,chCS] = max(abs(y));
      bh = mean(x(1:chCS));
      wh = mean(x(chCS+1:end)) - bh;
      xST = [x(1:chCS); x(chCS+1:end) - wh];
      for bs = 1:B
        bsInd = prmMat(:,randperm(N));
        xSTbsCUSUM = xST(bsInd);
        xSTbsCUSUM = xSTbsCUSUM(:);
        ySTbsCUSUM = cumsum(xSTbsCUSUM-mean(xSTbsCUSUM));
        SbsCUSUM(bs) = max(abs(ySTbsCUSUM));
      end
      
      % test for CUSUM's CP significance
      model.chCS   = chCS;
      model.chCSBS = [];
      if Scs > quantile(SbsCUSUM,1-.05)
        model.chCSBS = chCS;
      end
      
      % save model (includes both CUSUM and PARCS results)
      gSTR = [num2str(idivide(sd*10,int32(10))), 'p', num2str(mod(sd*10,10))];
      DataDirSTR = '../data/demo1_figure2/';
      if ~exist(DataDirSTR,'dir'), mkdir(DataDirSTR); end
      save([DataDirSTR , 'model_sd',gSTR,'_c',num2str(cp,'%02d'),...
        '_',num2str(mod(r,R),'%03d'),'.mat'],'model');
    end
  end
  fprintf('realization %d\n',r);
end

% (c) 2018 Hazem Toutounji, Dept. Theoretical Neuroscience, Central Institute of
% Mental Health, Medical Faculty Mannheim of Heidelberg University