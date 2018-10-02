%% Comparing CUSUM and PARCS for AMOC and MA(2) noise
% This code generates data for Figure 3, as
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
T = 100; % number of time steps
x = zeros(T,1); % time series
Cs = 20:10:80; C = numel(Cs); % ground truth CPs
Ws = .7:.1:1.3 ; W = numel(Ws); % CP weight
kappa = [.7 -.5 .4]/.7; p = numel(kappa);

% block-permutation bootstrap parameters
B = 10000; % number of bootstrap samples
k = 1    ; % block size
N = T/k  ; % number of blocks

% for block-permutation bootstrap of CUSUM
prmMat = reshape(1:T,k,N);
SbsCUSUM = zeros(B,1);

Err = zeros(T,1);

for r = 1:R
  for sd = [.7 1] % two levels of noise
    for c = 1:C
      for s = 1:W
        w = Ws(s);
        cp  = Cs(c);
        
        % generate MA(2) noise
        e = sd*randn(T+p-1,1);
        for t = 1:T
          Err(t) = kappa * e((1:p)+t-1);
        end
        
        % generate time series
        x(1:cp) = b + Err(1:cp);
        x(cp+1:end) = b + w + Err(cp+1:end);
        
        
        % run PARCS_1
        model = parcs(x,1);
        model.c = cp;
        model.b = b;
        model.w = w;
        model.kappa = kappa;
        model.sd = sd;
        
        % specify MA order (Algorithm 2) and block size k
        model.r = nan(10);
        y0 = model.y - model.yh(:,:,end);
        x0 = bsxfun(@plus,[y0(1,:);diff(y0)],mean(model.x));
        for k = 1:10
          xtdt = x0((k+1):T) - mean(x0((k+1):T));
          xt   = x0(1:(T-k)) - mean(x0(1:(T-k)));
          rho  = mean(xtdt.*xt)/(std(xtdt)*std(xt));
          prc  = normcdf(rho,-1/(T-k),sqrt(1/(T-k)));
          model.r(k) = rho;
          if prc >= .025 && prc <= .975
            break;
          end
        end
        model.r = model.r(~isnan(model.r));
        
        % run block-permutation bootstrap on PARCS_1
        model = bpb4parcs(model,k,B,.05);
        
        % CUSUM-transformed time series
        y = model.y;
        
        % specify MA order given CUSUM's CP and block size k
        for k = 1:10
          xtdt = x0((k+1):T) - mean(x0((k+1):T));
          xt   = x0(1:(T-k)) - mean(x0(1:(T-k)));
          rho  = mean(xtdt.*xt)/(std(xtdt)*std(xt));
          prc  = normcdf(rho,-1/(T-k),sqrt(1/(T-k)));
          
          if prc >= .025 && prc <= .975
            break;
          end
        end
        
        % bootstrap CUSUM
        [Scs,chCS] = max(abs(y));
        
        bh = mean(x(1:chCS));
        wh = mean(x(chCS+1:end)) - bh;
        x0 = [x(1:chCS); x(chCS+1:end) - wh];
        N = ceil(T/k);
        prmMat = reshape(1:N*k,k,N);
        for bs = 1:B
          bsInd = prmMat(:,randperm(N));
          bsInd = bsInd(bsInd<=T);
          xSTbsCUSUM = x0(bsInd);
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
        wSTR = [num2str(idivide(w*10,int32(10))), 'p', num2str(mod(w*10,10))];
        gSTR = [num2str(idivide(sd*10,int32(10))), 'p', num2str(mod(sd*10,10))];
        DataDirSTR = '../data/demo2_figure3/';
        if ~exist(DataDirSTR,'dir'), mkdir(DataDirSTR); end
        save([DataDirSTR , 'model_sd',gSTR,'_w',wSTR,'_c',num2str(cp,'%02d'),...
          '_',num2str(mod(r,R),'%03d'),'.mat'],'model');
        
      end
    end
  end
  fprintf('realization %d\n',r);
end

% (c) 2018 Hazem Toutounji, Dept. Theoretical Neuroscience, Central Institute of
% Mental Health, Medical Faculty Mannheim of Heidelberg University