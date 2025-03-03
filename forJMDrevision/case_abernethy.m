%% NOTES
% fixed importance sampling code on 08262016
% on 08262016, showed good gisa result when theta=100, see gisa_s10000_inq100_n1000_comp1_theta100_nt2455_08262016
% but theta=1 does not work, which is reasonable?

% theta = 100, s = 1e4, inq = 100, GISA DONE need expected values
% theta = 100, s = 1e5, inq = 100, GISA run once
% theta = 100, s = 1e4, inq = 10, GISA running
% theta = 1, s = 1e4, inq = 100, GISA 
% theta = 1, s = 1e4, Abernethy
% theta = 100, s = 1e4, Abernethy running

% parpool(4);
load('../basedata.mat');
addpath('../../Tools/liblinear');

c = Xf(:,26:30)*price'-cv; %price - cost
TEST = 20;
MAX_ITER = 1000;
prob_set = cell(TEST,1);
pairs_set = cell(TEST,1);
dx_set = cell(TEST,1);
partworths_set = cell(TEST,1);
target_best_set = zeros(TEST,1);
strategy_set = cell(TEST,1); % 1: determinant, 2: most probable 
expected_value_set = cell(TEST,1);
cond_set = cell(TEST,1);
s = 1e4;
inq = 100;
W0 = mvnrnd(zeros(1,d),eye(length(w)),s);
XID = 1:30; 
XID(5:5:30)=[];
sigma = 1e-36;
Dw = eye(30)*sigma; % randomness in user choices

nt = size(Xf,1); % number of testing object
% nt = 10;

theta = 100;
wtrue = w*theta;
wtrue(1:5) = wtrue(1:5)-wtrue(5);
wtrue(6:10) = wtrue(6:10)-wtrue(10);
wtrue(11:15) = wtrue(11:15)-wtrue(15);
wtrue(16:20) = wtrue(16:20)-wtrue(20);
wtrue(21:25) = wtrue(21:25)-wtrue(25);
wtrue(26:30) = wtrue(26:30)-wtrue(30);

num_competitor = 1;


for test = 1:TEST
    rng(test);
    fprintf('\n%%%%%%%% test number %d %%%%%%%%%%',test);
    
    % Calculate true distribution under sigma
    Wtrue = wtrue*theta;
    util = (Wtrue*Xf(1:nt,:)');
    competitors = randperm(nt,num_competitor);
    util_competitor = util(:,competitors);
    util_competitor_all = kron(util_competitor,ones(1,nt));
    util_all = repmat(util,1,num_competitor);
    exp_delta_util = exp(-util_all+util_competitor_all);
    exp_sum_delta_util = exp_delta_util*sparse(repmat(eye(nt),num_competitor,1))+1;
    obj_app = bsxfun(@plus,-log(exp_sum_delta_util),log(c'));
    
    best = bsxfun(@eq, obj_app, max(obj_app,[],2));
    best = best.*util;
    best(best==0) = -1e9;
    best = bsxfun(@eq, best, max(best,[],2));
    target_dist = sum(best/s,1)';
    [~,target_best] = max(target_dist,[],1);
    %     target_best = 1701;
    target_best_set(test) = target_best;

    queryID = dXID;
    nx = size(Xf,1);
    DX = dXID + dXID';
    probability_obj_set = zeros(nt,MAX_ITER);
    partworths = zeros(d,MAX_ITER);
    dxs = zeros(MAX_ITER,length(XID));
    conds = zeros(1,MAX_ITER);
    pairs = zeros(MAX_ITER,2);
    expected_values = zeros(nt,MAX_ITER);
    %     probability_obj_set = [prob_set{test},zeros(2455,500)];
    %     partworths = [partworths_set{test},zeros(30,500)];
    %     pairs = [pairs_set{test};zeros(500,2)];
    nq = 1;
    nquery = sum(sum(queryID>0));
    queryList = ones(nquery,1);
        
    % find the next query
    while  nq<=MAX_ITER
        % calculate A
        A = dX(DX(sub2ind(size(DX),pairs(1:nq-1,1),pairs(1:nq-1,2))),XID);
        A = bsxfun(@times,A,sign(pairs(1:nq-1,2)-pairs(1:nq-1,1)));
        [probability_obj,W,I,unique_I,w0,C,weights,expected_value] = ...
            appObjDistribution(s,length(XID),W0(:,XID),Xf(1:nt,XID),c(1:nt,:),A,competitors);
        
        util = Xf(1:nt,XID)*w0';
        obj_app = bsxfun(@plus,-log(1+exp(-util)),log(c(1:nt,:)));

        if nq>1
            tic
            As = bsxfun(@times, A, sqrt(exp(A*w0'))./(1+exp(A*w0')));
            As(isnan(As))=0;
            Sigma_inv = (eye(length(XID))/C+(As'*As));
            conds(nq) = cond(Sigma_inv);
            
            [~,guess] = max(probability_obj);
            [~,guess_expect] = max(expected_value);
            fprintf('iter: %d, truth: %d (%f), guess: %d, max value cand.: %d, corr: %0.2f, norm: %0.2f  \n',...
                nq, target_best, probability_obj(target_best), guess(1),...
                guess_expect(1), corr(w0',wtrue(XID)'), norm(w0'-wtrue(XID)'));
            
            B = (eye(size(w0,2))-w0'*w0/(w0*w0'))*Sigma_inv;
            [V,D] = eig(B);
            e = diag(D);
            v = V(:,e==min(e(e>1e-12)));
            v = v(:,1);
            s1 = dX(:,XID)*w0';
            s2 = (dX(:,XID)*v)./sqrt(sum(dX(:,XID).^2,2));

            unsampled = unique(queryID);
            unsampled(unsampled==0)=[];
%                 unsampled = setdiff(unsampled, zerorow);
            s2_id = find(abs(s2(unsampled))==max(abs(s2(unsampled))));
            s1_id = find(abs(s1(unsampled(s2_id)))==min(abs(s1(unsampled(s2_id)))));
            id = unsampled(s2_id(s1_id(1)));
            [a1,a2] = find(dXID==id);
            options = [a1,a2];
            toc
            aaa=1;
        else
            options = randperm(nx,2); % random initial query
        end

        if isempty(options)
            probability_obj_set(:,nq:end) = repmat(probability_obj,1,MAX_ITER-nq+1);
            expected_values(:,nq:end) = repmat(expected_value,1,MAX_ITER-nq+1);
            partworths(XID,nq:end) = repmat(w0',1,MAX_ITER-nq+1);
            nq = MAX_ITER+1;
        else
            options = options(1,:);
            queryID(min(options),max(options))=0;


            w_now = mvnrnd(wtrue,Dw,1)'*theta;
            u = (Xf(options(1),:)-Xf(options(2),:))*w_now;
            if rand()<1/(1+exp(-u))
                pairs(nq,:) = [options(1),options(2)];
            else
                pairs(nq,:) = [options(2),options(1)];
            end

            probability_obj_set(:,nq) = probability_obj;
            expected_values(:,nq) = expected_value;
            partworths(XID,nq) = w0';
            nq = nq+1;               
        end

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    prob_set{test} = probability_obj_set;
    pairs_set{test} = pairs;
    partworths_set{test} = partworths;
    cond_set{test} = conds;
    expected_value_set{test} = expected_values;
end
% save(['abernethy_s',num2str(s),'_n',num2str(MAX_ITER),...
%     '_comp',num2str(num_competitor),'_theta',num2str(theta),...
%     '_nt',num2str(nt),'_09042016.mat'],...
%     'pairs_set','prob_set','partworths_set','target_best_set','cond_set','expected_value_set','-v7.3');