%% Evaluating PARCS and block-size specification for 2 CPs and MA(2) noise
% This code generates data for Figure 6, as
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
T = 100; % number of time steps
t = (1:T)'; % time index
x = zeros(T,1); % time series
c = [20 60];
sd = .7; % s.d. of white Gaussian noise
kappa = [.7 -.5 .4]/sd; p = numel(kappa);

% block-permutation bootstrap parameters
B = 10000; % number of bootstrap samples

Err = zeros(T,1);

for r = 1:R
  for d = 1:3 % scenario
    
    % generate MA(2) noise
    e = sd*randn(T+p-1,1);
    for t = 1:T
      Err(t) = kappa * e((1:p)+t-1);
    end
    
    % generate time series
    w = W(d,:);
    x(1:c(1)) = b + Err(1:c(1));
    x(c(1)+1:c(2)) = b + w(1) + Err(c(1)+1:c(2));
    x(c(2)+1:end)  = b + w(1) + w(2) + Err(c(2)+1:end);
    
    % run PARCS_3, followed by block-permutation bootstrap
    model = parcs(x,3);
    model.c = c;
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
    
    model = bpb4parcs(model,k,B,.05);
    
    % save PARCS model
    DataDirSTR = '../data/demo5_figure6/';
    if ~exist(DataDirSTR,'dir'), mkdir(DataDirSTR); end
    save([DataDirSTR , 'model_scenario',num2str(d),...
      '_',num2str(mod(r,R),'%03d'),'.mat'],'model');
    
  end
  fprintf('realization %d\n',r);
end

% (c) 2018 Hazem Toutounji, Dept. Theoretical Neuroscience, Central Institute of
% Mental Health, Medical Faculty Mannheim of Heidelberg University