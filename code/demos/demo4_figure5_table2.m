%% Comparing CUSUM and PARCS for 2 CPs and white Gaussian noise
% This code generates data for Figure 5 and Table 2, as
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
W = [1 2; 2 -1; 2 1];
C = [20 60];
sig = 1; % s.d. of white Gaussian noise

% block-permutation bootstrap parameters
B = 10000; % number of bootstrap samples
k = 1    ; % block size

% for block-permutation bootstrap of CUSUM
SbsCUSUM = zeros(B,1);

for r = 1:R
  for T = [26 , 50 , 100] % number of time steps
    c = round(C*T/100); % change points
    N = T/k  ; % number of bloks
    prmMat = reshape(1:T,k,N);
    for d = 1:3 % scenario
      
      % generate time series
      w = W(d,:); x = zeros(T,1);
      x(1:c(1)) = b + sig*randn(c(1),1);
      x(c(1)+1:c(2)) = b + w(1) + sig*randn(c(2)-c(1),1);
      x(c(2)+1:end)  = b + w(1) + w(2) + sig*randn(T-c(2),1);
      
      % run PARCS_3, followed by block-permutation bootstrap
      model = parcs(x,3);
      model.c = c;
      model.b = b;
      model.w = w;
      model.sd = sig;
      model = bpb4parcs(model,k,B,.30);
      
      % -------------------------------------------
      % CUSUM with binary segmentation (terminated
      % after at most one partitioning; 3CPs max)
      % -------------------------------------------
      
      % CUSUM trnasformation on full series
      x1 = model.x;
      y1 = model.y;
      
      % bootstrap CUSUM on full series
      model.chBinSeg   = nan(3,1);
      model.chBinSegBS = nan(3,1);
      [Scs,chCS] = max(abs(y1));
      prmMat = reshape(1:T,k,N);
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
      model.chBinSeg(1) = chCS;
      if Scs > quantile(SbsCUSUM,1-.05)
        model.chBinSegBS(1) = chCS;
      end
      
      % define candidate c_21 and c_22 for binary segmentation
      x2 = x1(1:model.chBinSeg(1)-1);
      x3 = x1(model.chBinSeg(1)+1:end);
      T2 = numel(x2); N2 = T2/k;
      T3 = numel(x3); N3 = T3/k;
      y2 = cumsum(x2-mean(x2));
      y3 = cumsum(x3-mean(x3));
      [~,model.chBinSeg(2)] = max(abs(y2));
      [~,model.chBinSeg(3)] = max(abs(y3));
      model.chBinSeg(3) = model.chBinSeg(3) + model.chBinSeg(1);
      
      % if c_1 confirmed, bootstrap c_21 and c_22
      if ~isnan(model.chBinSeg(1))
        
        % bootstrap for c_21
        [Scs,chCS] = max(abs(y2));
        prmMat = reshape(1:T2,k,N2);
        x = x2;
        bh = mean(x(1:chCS));
        wh = mean(x(chCS+1:end)) - bh;
        xST = [x(1:chCS); x(chCS+1:end) - wh];
        for bs = 1:B
          bsInd = prmMat(:,randperm(N2));
          xSTbsCUSUM = xST(bsInd);
          xSTbsCUSUM = xSTbsCUSUM(:);
          ySTbsCUSUM = cumsum(xSTbsCUSUM-mean(xSTbsCUSUM));
          SbsCUSUM(bs) = max(abs(ySTbsCUSUM));
        end
        if Scs > quantile(SbsCUSUM,1-.05)
          model.chBinSegBS(2) = model.chBinSeg(2);
        end
        
        % bootstrap for c_22
        [Scs,chCS] = max(abs(y3));        
        prmMat = reshape(1:T3,k,N3);
        x = x3;
        bh = mean(x(1:chCS));
        wh = mean(x(chCS+1:end)) - bh;
        xST = [x(1:chCS); x(chCS+1:end) - wh];
        for bs = 1:B
          bsInd = prmMat(:,randperm(N3));
          xSTbsCUSUM = xST(bsInd);
          xSTbsCUSUM = xSTbsCUSUM(:);
          ySTbsCUSUM = cumsum(xSTbsCUSUM-mean(xSTbsCUSUM));
          SbsCUSUM(bs) = max(abs(ySTbsCUSUM));
        end
        if Scs > quantile(SbsCUSUM,1-.05)
          model.chBinSegBS(3) = model.chBinSeg(3);
        end
      end
      
      % save model (includes both binary segmentation and PARCS results)
      DataDirSTR = '../data/demo4_figure5_table2/';
      if ~exist(DataDirSTR,'dir'), mkdir(DataDirSTR); end
      save([DataDirSTR , 'model_scenario',num2str(d),'_T',num2str(T,'%03d'),...
        '_',num2str(mod(r,R),'%03d'),'.mat'],'model');
      
    end
  end
  fprintf('realization %d\n',r);
end

% (c) 2018 Hazem Toutounji, Dept. Theoretical Neuroscience, Central Institute of
% Mental Health, Medical Faculty Mannheim of Heidelberg University