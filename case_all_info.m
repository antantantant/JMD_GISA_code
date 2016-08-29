% parpool('local',12);
% % clear;
% % 
% % test case: all comparisons known
% load basedata.mat;
addpath('C:\Users\Max Yi Ren\Documents\Research\Code\Tools\liblinear\matlab');
target = load('target_distribution.mat','probability_obj');
target = target.probability_obj;
target = target(:);

TEST = 1;

c = Xf(:,26:30)*price'-cv;
prob_set = cell(TEST,1);
pairs_set = cell(TEST,1);
s = 1e5;
inq = 10;
alg = 'profit';

for test = 1:TEST
    rng(test);
    fprintf('\n%%%%%%%% test number %d %%%%%%%%%%',test);
    [prob,pairs] = ...
        appGGBS(s,d,Xf,dX,dXID,c,inq,[],alg,dXID,[],w',target,1);
%     [prob,pairs] = bestProb(s,d,Xf,dX,dXID,c,inq,[],alg,dXID,[],w',target,1);

    prob_set{test} = prob;
    pairs_set{test} = pairs;
end
save(['profit_s',num2str(s),'_inq',num2str(inq),...
    '.mat'],'prob_set','pairs_set');

% save(['bestprob_s',num2str(s),'_inq',num2str(inq),...
%     '.mat'],'prob_set','pairs_set');