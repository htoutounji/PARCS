function model = bpb4parcs(model,k,B,alpha)
% Block-Permutation Bootstrap for PARCS, Algorithm 3, as
% described in Toutounji and Durstewitz (2018) Front Neuroinform.
% Please refer to 10.3389/fninf.2018.00067 for details
% INPUT:
%  model: PARCS model structure;
%  k    : block size for bootstrap
%  B    : number of bootstrap samples
%  alpha: nominal significance level
% OUTPUT:
%  model: PARCS model structure; additional fields
%    model.BBS: regression coefficients for bootstrapped PARCS model
%    model.yhBS: bootstrapped PARCS model
%    model.chBS: bootstrapped, ranked CPs
%    model.HBS: splines in the bootstrapped PARCS model
%    model.tst: test statistic for CP_(1:M)

M = model.M;
N = model.N;
T = model.T;

yh = reshape(model.yh,[T,N,M+1]);
model.k    = k;
model.BBS  = reshape(model.B,[2*M+1,N,M+1]) ;
model.yhBS = yh;

model.chBS = model.ch;
model.HBS  = model.H ;

cconf = true(2*M,1);

% computing H0-comform time series
y0 = model.y - yh(:,:,end);
x0 = bsxfun(@plus,[y0(1,:);diff(y0)],mean(model.x));

model.tst = zeros(M,1);
tstBS = zeros(B,1);
Nb   = ceil(T/k);
prmMat = reshape(1:k*Nb,k,Nb);

for m = 1:M

  % regressing out significant CPs
  ind = 1:2*(m-1)+1;
  ind = ind([true;cconf(ind(1:end-1))]);
  ydfh = model.H(:,ind)*((model.H(:,ind)'*model.H(:,ind))\(model.H(:,ind)'*model.y));
  y1 = model.y - ydfh;

  % computing test statistic for CP_m
  ind = [1,2*m:2*M+1];
  B1 = (model.H(:,ind)'*model.H(:,ind))\(model.H(:,ind)'*y1);
  model.tst(m) = mean(abs(B1(2,:)+B1(3,:)));
  
  % estimating EDF through block-permutation bootstrap
  for b = 1:B
    bsInd = prmMat(:,randperm(Nb));
    bsInd = bsInd(bsInd<=T);
    
    x0bs = x0(bsInd,:);
    y0bs = cumsum(bsxfun(@minus,x0bs,mean(x0bs)));
    
    B0 = (model.H(:,ind)'*model.H(:,ind))\(model.H(:,ind)'*y0bs);

    tstBS(b) = mean(abs(B0(2,:)+B0(3,:)));
  end
  
  % test for significance
  if model.tst(m) <= quantile(tstBS,1-alpha)
    cconf(2*m-1:2*m) = false;
    model.chBS(m) = nan;
    model.HBS(:,2*m:2*m+1) = nan;
  end
end

% remove nan values (insignificant CPs and their splines)
c = 1;
for m = 1:M
  if isnan(model.HBS(1,2*c))
    model.HBS = [model.HBS(:,1:2*(c-1)+1) , model.HBS(:,2*(c+1):end)];
    model.chBS(c) = [];
  else
    c = c+1;
  end
end
    
% compute regression coefficients and final PARCS model
model.BBS  = (model.HBS'*model.HBS)\(model.HBS'*model.y);
model.yhBS = model.HBS*model.BBS;

% for univariate time series
if N == 1
  model.B = squeeze(model.B);
  model.yh = squeeze(model.yh);
end

% (c) 2018 Hazem Toutounji, Dept. Theoretical Neuroscience, Central Institute of
% Mental Health, Medical Faculty Mannheim of Heidelberg University