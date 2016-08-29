% parpool('local',12);
% load('..\basedata.mat');
addpath('..\..\Tools\liblinear\matlab');

c = Xf(:,26:30)*price'-cv; %price - cost
TEST = 1;
MAX_ITER = 500;
prob_set = cell(TEST,1);
pairs_set = cell(TEST,1);
partworths_set = cell(TEST,1);
s = 1e4;
inq = 10;
W0 = mvnrnd(zeros(1,d),eye(length(w)),s);
XID = 1:30; 
XID(5:5:30)=[];
sigma = 1e-6;
Dw = eye(30)*sigma; % randomness in user choices

% nt = size(Xf,1); % number of testing object
nt = 100;

% Calculate true distribution under sigma
Wtrue = mvnrnd(w,Dw,s);
theta = 1;
Wtrue = theta*Wtrue;
util = (Wtrue*Xf(1:nt,:)');
obj_app = bsxfun(@plus,-log(1+exp(-util)),log(c(1:nt,:)'));
best = bsxfun(@eq, obj_app, max(obj_app,[],2));
best = best.*util;
best(best==0) = -1e9;
best = bsxfun(@eq, best, max(best,[],2));
target_dist = sum(best/s)';
[~,target_best] = max(target_dist);

for test = 1:TEST
    rng(test);
    fprintf('\n%%%%%%%% test number %d %%%%%%%%%%',test);
    
    wtrue = theta*w';
    queryID = dXID;
    nx = size(Xf,1);
    DX = dXID + dXID';
    probability_obj_set = zeros(nt,MAX_ITER);
    nq = 100000;
    dX_shuffle = dX(randperm(size(dX,1)),:);
    
    dX_t = dX_shuffle(1:nq,:);
    y = rand(nq,1)<(1./(1+exp(-dX_t*wtrue)));
    A = bsxfun(@times,dX_t,y);


    [probability_obj,W,I,unique_I,w0,C,weights] = ...
        appObjDistribution(s,length(XID),W0(:,XID),Xf(1:nt,XID),c(1:nt,:),A(:,XID));
        
    util = Xf(1:nt,XID)*w0';
    obj_app = bsxfun(@plus,-log(1+exp(-util)),log(c(1:nt,:)));
        
        
    [~,guess] = max(probability_obj);
    fprintf('iter: %d, truth: %d (%f), guess: %d, point guess: %d \n',...
        nq, target_best, probability_obj(target_best), guess, find(obj_app==max(obj_app),1));
        
end
