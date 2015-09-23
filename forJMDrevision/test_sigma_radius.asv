%%% This code tests how large information matrix should be to get the
%%% correct answer

% matlabpool open;
load('..\basedata.mat');
addpath('..\..\Research\Code\Tools\liblinear\matlab');
rng(1);
c = Xf(:,26:30)*price'-cv; %price - cost
s = 1e4;
inq = 100;
W0 = mvnrnd(zeros(1,d),eye(length(w)),s);
XID = 1:30; 
XID(5:5:30)=[];
sigma = 1e-36;
Dw = eye(30)*sigma; % randomness in user choices

% nt = size(Xf,1); % number of testing object
nt = 10;

theta = 1;
wtrue = w*theta;
wtrue(1:5) = wtrue(1:5)-wtrue(5);
wtrue(6:10) = wtrue(6:10)-wtrue(10);
wtrue(11:15) = wtrue(11:15)-wtrue(15);
wtrue(16:20) = wtrue(16:20)-wtrue(20);
wtrue(21:25) = wtrue(21:25)-wtrue(25);
wtrue(26:30) = wtrue(26:30)-wtrue(30);

num_competitor = 1;

% Calculate true distribution under sigma
Wtrue = wtrue*theta;
util = (Wtrue*Xf(1:nt,:)');
competitors = randperm(nt,num_competitor);
util_competitor = util(:,competitors);
util_competitor_all = kron(util_competitor,ones(1,nt));
util_all = repmat(util,1,num_competitor);
exp_delta_util = exp(-util_all+util_competitor_all);
exp_sum_delta_util = exp_delta_util*repmat(eye(nt),num_competitor,1)+1;
obj_app = bsxfun(@plus,-log(exp_sum_delta_util),log(c(1:nt,:)'));
best = bsxfun(@eq, obj_app, max(obj_app,[],2));
best = best.*util;
best(best==0) = -1e9;
best = bsxfun(@eq, best, max(best,[],2));
target_dist = sum(best/s,1)';
[~,target_best] = max(target_dist,[],1);
%     target_best = 1701;
target_best_set = target_best;

    
queryID = dXID;
nx = size(Xf,1);
DX = dXID + dXID';


X = Xf(1:nt,XID);
A = dxs(1:nq-1,:);
[w0, C, weights] = crossvalidation(sparse(A),ones(size(A,1),1)); 
As = bsxfun(@times, A, sqrt(exp(A*w0'))./(1+exp(A*w0')));
As(isnan(As))=0;
Sigma_inv = (eye(length(XID))/C+(As'*As));
W = mvnrnd(wtrue(XID)/norm(wtrue)*norm(w0),0.1*inv(Sigma_inv),s);
util = (W*X');
num_competitior = length(competitors);
util_competitor = util(:,competitors);
util_competitor_all = kron(util_competitor,ones(1,size(X,1)));
util_all = repmat(util,1,num_competitior);
exp_delta_util = exp(-util_all+util_competitor_all);
exp_sum_delta_util = exp_delta_util*repmat(eye(size(X,1)),num_competitior,1)+1;
obj_app = bsxfun(@plus,-log(exp_sum_delta_util),log(c(1:nt,:)'));

best = bsxfun(@eq, obj_app, max(obj_app,[],2));
best = best.*util;
best(best==0) = -1e9;
best = bsxfun(@eq, best, max(best,[],2));
f = (kron(weights,ones(1,size(W0,1)))*best/s)';

% each line of best needs to have only one 1
if any(sum(best,2)>1)
    error = 1;
end

I = best*(1:size(best,2))';
%     obj = bsxfun(@times,1./(1+exp(-W*X(keep,:)')),c(keep)');
%     f(keep) = sum(bsxfun(@eq, obj, max(obj,[],2))/s)';

unique_I = unique(I);
unique_I(unique_I>length(c))=[];

probability_obj = f;
[~,guess] = max(probability_obj);
fprintf('truth: %d (%f), guess: %d\n',...
    target_best, probability_obj(target_best), guess);
