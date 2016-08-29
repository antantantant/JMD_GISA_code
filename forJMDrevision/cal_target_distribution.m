load('..\basedata.mat');
addpath('C:\Users\Max Yi Ren\Documents\Research\Code\Tools\liblinear\matlab');

% %% preprocessing
[n,p] = size(dX);
dXunique = unique(dX, 'rows'); %keep unique rows
A = bsxfun(@times, dXunique, (abs(dXunique*w')>1e-6).*sign(dXunique*w')); %flip the pair if the second has higher utility
A(sum(A.^2,2)<1e-6,:) = []; %remove all-zero rows

%% get w distribution
% [w,C] = crossvalidation(sparse([A;-A]), [ones(size(A,1),1);-ones(size(A,1),1)]);
% model = train([ones(size(A,1),1);-ones(size(A,1),1)],...
%                     sparse([A;-A]),'-s 0 -c 1e6 -e 1e-6 -q');
% w_all = model.w;

load fullinfo_w.mat; %get w and C

%% sample according to the distribution of w

sample_size = 1e4;
feature_id = 1:30;
C = 1;
% use metropolis-hasting to draw from the posterior distribution of w
logpdf_h = @(x) -sum(log(1+exp(-A*x'))) - 1/2/C*x*x';
proprnd_h = @(x) x + (rand(1,p) - 0.5);

model = train(ones(size(A,1),1),...
        sparse(A), sprintf('-s 0 -e 1e-6 -q -c %f', C));
w0 = model.w;
tic
W1 = mhsample(w0,sample_size,'logpdf',logpdf_h,'proprnd',proprnd_h,...
    'symmetric',true,'burnin',sample_size);
toc

% %use codegen to convert this code into mex
% start_t = coder.typeof(0.1,[1, Inf],[0 1]);
% nsamples_t = coder.typeof(10);
% burnin_t = coder.typeof(10);
% C_t = coder.typeof(0.1);
% A_t = coder.typeof(0.1,[Inf, Inf],[1 1]);
% codegen myMetropolisHasting -args ...
%     {start_t, nsamples_t, burnin_t, A_t, C_t};
tic
W2 = myMetropolisHasting_mex(zeros(1,p),sample_size*2,sample_size,A(1:50,:),C);
toc


% % calculate probability density at all samples
% pw = exp(-C*sum(log(1+exp(-A*W'))) - 1/2*sum(W'.^2));
% pw = pw/sum(pw);


c = Xf(:,26:30)*price'-cv;
util = (W1*Xf(:,feature_id)');
obj_app = bsxfun(@plus,-log(1+exp(-util)),log(c'));
best = bsxfun(@eq, obj_app, max(obj_app,[],2));
best = best.*util;
best = bsxfun(@eq, best, max(best,[],2));
f_app = sum(best/sample_size)';
find(f_app>0)





CC = 2e-5;
Scale = 1;
addpath('..\..\Tools\liblinear\matlab');
% feature_id = [2:5,7:10,12:15,17:20,22:25,27:30];
feature_id = 1:30;
As = A(:,feature_id);

model = train(ones(size(As,1),1),...
        sparse(As), sprintf('-s 0 -e 1e-6 -q -c %f', CC));
w0 = model.w;
Ad = bsxfun(@times, As, sqrt(exp(As*w0'))./(1+exp(As*w0')));
Sigma = (eye(length(w0))/CC+(Ad'*Ad))*Scale;

% % Toubia method
% w0 = (As'*As+1/CC*eye(size(As,2)))\As'*ones(size(As,1),1);
% Sigma = As'*As+1/CC*eye(size(As,2));

sample_size = 1e4;
W_app = mvnrnd(w0,inv(Sigma),sample_size);


util = (W_app*Xf(:,feature_id)');
obj_app = bsxfun(@plus,-log(1+exp(-util)),log(c'));
best = bsxfun(@eq, obj_app, max(obj_app,[],2));
best = best.*util;
best = bsxfun(@eq, best, max(best,[],2));
f_app = sum(best/sample_size)';
find(f_app>0)


%% consider a distribution of the true partworth
true_var = 0.01;
W_true = mvnrnd(w, true_var*eye(length(w)), sample_size);
util = (W_true*Xf(:,feature_id)');
obj_app = bsxfun(@plus,-log(1+exp(-util)),log(c'));
best = bsxfun(@eq, obj_app, max(obj_app,[],2));
best = best.*util;
best = bsxfun(@eq, best, max(best,[],2));
f_app = sum(best/sample_size)';
find(f_app>0)