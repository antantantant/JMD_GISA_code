% main entrance for detc2014b on GGBS for nonlinear objective

addpath('..\..\..\Code\Tools\liblinear\matlab');
% set random seed
rng(0);

%%%%%%%%%%%%%%%%%%%%% 2d test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx = 10; % number of products
d = 10; % number of product attributes
X = (rand(nx,d)-0.5)*2; % product attributes
c = rand(nx,1); % price - cost
w = (rand(d,1)-0.5)*2; % true w

% find true best product
f = 1./(1+exp(-X*w)).*c;
true_best = find(f==max(f));

nq = nx*(nx-1)/2; % number of queries
dX = zeros(nq,d); % attribute difference data
dXID = zeros(nx,nx); % pair to query id mapping
count = 1;
for i = 1:nx
    for j = i+1:nx
        dXID(i,j) = count;
        dX(count,:) = X(i,:) - X(j,:);
        count = count + 1;
    end
end

% experiment
par.X = X;
par.c = c;
par.dX = dX;
par.dXID = dXID;
par.queryID = dXID;
par.nx = nx;
par.d = d;
par.nq = nq;
par.candidate = 1:par.nx;

% parameters for sampling
par.s = 1e5;
par.W = (rand(par.s,par.d)-0.5)*2;
pref = par.W*par.X';
area = (sum(bsxfun(@eq, pref, max(pref,[],2)))/par.s)';
const = 1/sum(area>0)./area;
const(area==0) = 0;
par.h = bsxfun(@eq, pref, max(pref,[],2))*const;

% active learning
[prob,pairs] = appGGBS(par,[],w,1);

% find best product
best = find(prob==max(prob));