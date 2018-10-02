function model = parcs(x,M)
% PARCS: Paired Adaptive Regressors for Cumulative Sum, Algorithm 1, as
% described in Toutounji and Durstewitz (2018) Front Neuroinform.
% Please refer to 10.3389/fninf.2018.00067 for details
% INPUT:
%  x: time series matrix of size TxN
%  M: number of change points
% OUTPUT:
%  model: PARCS model structure:
%    model.x: x
%    model.y: CUSUM-transformed time series
%    model.L: forward upper bound order
%    model.M: final upper-bound model complexity (following pruning)
%    model.N: number of covariates
%    model.T: number of time steps
%    model.hp: right-side spline dictionary
%    model.hm: left-side spline dictionary
%    model.ch: ranked CPs
%    model.mse: mean-square-error for PARCS_(0:M)
%    model.B: regression coefficients for PARCS_(0:M) models (size: 2M+1,N,M+1)
%    model.H: splines in the PARCS_M model
%    model.yh: PARCS_(0:M) models

warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');

% compute CUSUM transformation
y = cumsum(bsxfun(@minus,x,mean(x)));

% length T, dimensionality N and forward upper bound L
[T,N] = size(x); L = 3*M;

% initialising model
model.x = x;
model.y = y;
model.L = L;
model.M = M;
model.N = N;
model.T = T;

% generating splines library
Hpm = toeplitz(linspace(0,1,T));
hp  = tril(Hpm); model.hp = hp;
hm  = triu(Hpm); model.hm = hm;


cp  = nan(L,1); % change points
mse = nan(T,L);

%% forward stage

H = ones(T,2*L+1); % forward model splines
for m = 1:L
  for t = 2:T-1 % this loop can be optimised
    if ~ismember(t,cp(1:m-1))
      H(:,2*m)   = hp(:,t);
      H(:,2*m+1) = hm(:,t);
      
      B  = (H(:,1:2*m+1)' * H(:,1:2*m+1)) \ (H(:,1:2*m+1)' * y); % lse
      yh = H(:,1:2*m+1)*B;
      mse(t,m) = mean(mean((y-yh).^2));
    end
  end
  [~,cp(m)] = min(mse(:,m));
  H(:,2*m)   = hp(:,cp(m));
  H(:,2*m+1) = hm(:,cp(m));
end



%% backward stage

R = L - M; % number of removed CPs
for r = 1:R
  err = inf;
  for m = 1:L-r+1
    ind = setdiff(1:2*L-2*r+3,[2*m,2*m+1]);
    B = (H(:,ind)' * H(:,ind)) \ (H(:,ind)' * y);
    yh = H(:,ind)*B;
    mseR = mean((y(:)-yh(:)).^2);
    if mseR < err
      err     = mseR;
      remove  = ind;
      Cred = setdiff(cp,cp(m),'stable');
    end
  end
  H = H(:,remove);
  cp = Cred;
end
B = (H' * H) \ (H' * y);
yh = H*B;


%% ranking stage

model.ch  = zeros(M,1)      ;
model.mse = zeros(M+1,1)    ;
model.B   = nan(1+2*M,N,M+1);
model.H   = ones(T,2*M+1)   ;
model.yh  = zeros(T,N,M+1)  ;

model.mse(end) = mean((y(:)-yh(:)).^2);
for r = 1:M
  err = inf;
  for m = 1:M-r+1
    ind = setdiff(1:2*M-2*r+3,[2*m,2*m+1]);
    B = (H(:,ind)' * H(:,ind)) \ (H(:,ind)' * y);
    yh = H(:,ind)*B;
    mseR = mean((y(:)-yh(:)).^2);
    if mseR < err
      err     = mseR;
      remove  = ind;
      Cred = setdiff(cp,cp(m),'stable');
      Crem = cp(m);
    end
  end
  model.mse(end-r) = err;
  H = H(:,remove);
  cp = Cred;
  model.ch(end-r+1) = Crem;
end

%% Infer nested models

model.H(:,2:2:end) = hp(:,model.ch);
model.H(:,3:2:end) = hm(:,model.ch);

for m = 0:M
  model.B(1:2*m+1,:,m+1) = (model.H(:,1:2*m+1)'*model.H(:,1:2*m+1))\(model.H(:,1:2*m+1)'*y);
  model.yh(:,:,m+1) = model.H(:,1:2*m+1)*model.B(1:2*m+1,:,m+1);
end

% (c) 2018 Hazem Toutounji, Dept. Theoretical Neuroscience, Central Institute of
% Mental Health, Medical Faculty Mannheim of Heidelberg University