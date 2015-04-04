% parpool('local',12);
% % clear;
% % 
% % test case: all comparisons known
load('..\basedata.mat');
addpath('..\..\Tools\liblinear\matlab');
% target = load('..\target_distribution.mat','probability_obj');
% target = target.probability_obj;
% target = target(:);

target = sparse(zeros(size(Xf,1),1));
target(1704) = 1;

TEST = 1;

c = Xf(:,26:30)*price'-cv;
prob_set = cell(TEST,1);
pairs_set = cell(TEST,1);
partworths_set = cell(TEST,1);
s = 1e5;
W0 = mvnrnd(w,eye(length(w)),s);
alg = 'pref_only';
C = 1;

Dw = eye(30); % randomness in user choices

for test = 1:TEST
    rng(test);
    fprintf('\n%%%%%%%% test number %d %%%%%%%%%%',test);

%     [prob,pairs] = ...
%         GBS(s,d,Xf,dX,dXID,c,inq,[],alg,dXID,[],w',target,1);

%     half the version space, computationally cheaper version
    [prob,pairs,partworths] = ...
        halfV(s,d,Xf,dX,dXID,c,alg,dXID,W0,w',Dw,target,C);
    prob_set{test} = prob;
    pairs_set{test} = pairs;
    partworths_set{test} = partworths;
end

% save(['gbs_s',num2str(s),'_inq',num2str(inq),...
%     '.mat'],'prob_set','pairs_set');

% save(['toubia_s',num2str(s),'_C1_n1000_0220.mat'],...
%     'prob_set','pairs_set','partworths_set','-v7.3');