% parpool('local',12);
% % clear;
% % 
% % test case: all comparisons known
% load('..\basedata.mat');
addpath('..\..\Tools\liblinear\matlab');
% target = load('..\target_distribution.mat','probability_obj');
% target = target.probability_obj;
% target = target(:);

target = sparse(zeros(size(Xf,1),1));
target(1704) = 1;

TEST = 12;

c = Xf(:,26:30)*price'-cv;
% prob_set = cell(TEST,1);
% pairs_set = cell(TEST,1);
% partworths_set = cell(TEST,1);
s = 1e5;
inq = 1;
alg = 'profit';
true_var = 0; %variance of the true distribution
W0 = mvnrnd(w,eye(length(w)),s);
C = 1;

parfor test = 2:TEST
    rng(test);
    fprintf('\n%%%%%%%% test number %d %%%%%%%%%%',test);
    [prob,pairs,partworths] = ...
        appGGBS_PLUS(s,d,Xf,dX,dXID,c,inq,alg,dXID,W0,w',target,C, ...
        prob_set{test}, pairs_set{test}, partworths_set{test});
%     [prob,pairs] = bestProb(s,d,Xf,dX,dXID,c,inq,[],alg,dXID,[],w',target,1);

    prob_set{test} = prob;
    pairs_set{test} = pairs;
    partworths_set{test} = partworths;
end
% save(['profit_s',num2str(s),'_inq',num2str(inq),...
%     '_C1_n1000_bestonly_0228.mat'],'prob_set','pairs_set','partworths_set','-v7.3');

% save(['bestprob_s',num2str(s),'_inq',num2str(inq),...
%     '.mat'],'prob_set','pairs_set');