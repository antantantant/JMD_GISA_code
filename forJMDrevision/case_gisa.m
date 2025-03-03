%% NOTES
% fixed importance sampling code on 08262016
% on 08262016, showed good gisa result when theta=100, see gisa_s10000_inq100_n1000_comp1_theta100_nt2455_08262016
% but theta=1 does not work, which is reasonable?
% on 09042016, I found that in sampling, pw and qw are smaller than 1e-36
% in some cases, so removed these artificial terms, testing gisa,
% theta=100, s=1e4, inq=100 to see the effect of this.

% updated on 09052016
% theta = 100, s = 1e4, inq = 100, GISA DONE 09042016
% theta = 100, s = 1e3, inq = 100, GISA running on my office computer
% theta = 100, s = 1e5, inq = 100, GISA skip
% theta = 100, s = 1e4, inq = 10, GISA running on doi3
% theta = 1, s = 1e4, inq = 100, GISA running on doi3
% theta = 1, s = 1e4, inq = 100, Abernethy 
% theta = 100, s = 1e4, inq = 100, Abernethy running on my office computer

%%


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
CC_set = cell(TEST,1); 
expected_value_set = cell(TEST,1);
cond_set = cell(TEST,1);
s = 1e3;
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
    obj_app = bsxfun(@plus,-log(exp_sum_delta_util),log(c(1:nt,:)'));
    
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
    pairs = zeros(MAX_ITER,2);
    conds = zeros(1,MAX_ITER);
    expected_values = zeros(nt,MAX_ITER);
    minCC = zeros(MAX_ITER,1);
    strategy = zeros(MAX_ITER,1);
% %     probability_obj_set = prob_set{test};
%     partworths = partworths_set{test};
%     pairs = pairs_set{test};
    nq = 1;
    nquery = sum(sum(queryID>0));
        
    % find the next query
    while  nq<=MAX_ITER
        tic
        % calculate A
        A = dX(DX(sub2ind(size(DX),pairs(1:nq-1,1),pairs(1:nq-1,2))),XID);
        A = bsxfun(@times,A,sign(pairs(1:nq-1,2)-pairs(1:nq-1,1)));
        [probability_obj,W,I,unique_I,w0,C,weights,expected_value] = ...
            appObjDistribution(s,length(XID),W0(:,XID),Xf(1:nt,XID),c(1:nt,:),A,competitors);
        
        util = Xf(1:nt,XID)*w0';
        obj_app = bsxfun(@plus,-log(1+exp(-util)),log(c(1:nt,:)));
        
        if nq>1
            As = bsxfun(@times, A, sqrt(exp(A*w0'))./(1+exp(A*w0')));
            As(isnan(As))=0;
            Sigma_inv = (eye(size(A,2))/C+(As'*As));
            B = (eye(size(A,2))-(w0'*w0)/(w0*w0'))*(eye(size(A,2))/C+(As'*As));
            conds(nq) = cond(Sigma_inv);
            
            [~,guess] = max(probability_obj);
            [~,guess_expect] = max(expected_value);
            fprintf('iter: %d, truth: %d (%f), guess: %d, max value cand.: %d, corr: %0.2f, norm: %0.2f  \n',...
                nq, target_best, probability_obj(target_best), guess(1),...
                guess_expect(1), corr(w0',wtrue(XID)'), norm(w0'-wtrue(XID)'));
            
            [~,sort_obj] = sort(probability_obj,'descend');
            sort_obj_nonzero = sort_obj(1:sum(probability_obj>0));
            a = [];
            aid1 = [];
            count = 1;
            
            for i = 1:(length(sort_obj_nonzero)-1)
                for j = (i+1):length(sort_obj_nonzero)
                    if count <= inq/2
                        ii = min([sort_obj_nonzero(i),sort_obj_nonzero(j)]);
                        II = max([sort_obj_nonzero(i),sort_obj_nonzero(j)]);
                        if queryID(ii,II)>0
                            a = [a;[ii,II]];
                            aid1 = [aid1;queryID(ii,II)];
                            count = count + 1;
                        end
                    end
                end
            end
            
            if count <= inq/2
                sort_obj_zero = setdiff((1:nx)',sort_obj_nonzero);
                sort_obj_zero = sort_obj_zero(randperm(length(sort_obj_zero)));
                for i = 1:length(sort_obj_nonzero)
                    for j = 1:length(sort_obj_zero)
                        if count <= inq/2
                            ii = min([sort_obj_nonzero(i),sort_obj_zero(j)]);
                            II = max([sort_obj_nonzero(i),sort_obj_zero(j)]);
                            if queryID(ii,II)>0
                                a = [a;[ii,II]];
                                aid1 = [aid1;queryID(ii,II)];
                                count = count + 1;
                            end
                        end
                    end
                end
            end
%             aid = aid1; % test the query strategy with the best probablity
            
            unsampled = unique(queryID);
            unsampled(unsampled==0)=[];
            [V,D] = eig(B);
            e = diag(D);
            v = V(:,e==min(e(e>1e-12)));
            v = v(:,1);
            s2 = abs((dX(unsampled,XID)*v)./sqrt(sum(dX(unsampled,XID).^2,2)));
%             aid = unsampled(s2==max(s2));
%             aid = aid(randperm(length(aid),min([length(aid),inq])));
            [~,sort_id] = sort(s2,'descend');
            aid = unsampled(sort_id(1:inq/2));
            aid = [aid;aid1];
            
            actual_inq = length(aid);
            if actual_inq>0
                temp_set = zeros(nt,actual_inq);
                Q_remain = dX(aid,XID);
                u = W*Q_remain'>0; % get all utility signs
                u = bsxfun(@times, u, I); % label all positive signs with query numbers
                for j = 1:length(unique_I)
                    temp_set(unique_I(j),:) = (weights'*(u==unique_I(j)))/...
                        (weights'*(I==unique_I(j))+1e-99);
                end

                probability_query_pos = sum(bsxfun(@times,temp_set,probability_obj),1)';% probability of new query being positive

                rho_k = zeros(nt,actual_inq);
                rho = zeros(actual_inq,1);
                probability_query_pos(probability_query_pos==0)=1e-99;

                rho_k(temp_set>0.5) = temp_set(temp_set>0.5);
                rho_k(temp_set<=0.5) = 1-temp_set(temp_set<=0.5);
                rho(probability_query_pos>0.5) = probability_query_pos(probability_query_pos>0.5);
                rho(probability_query_pos<=0.5) = 1-probability_query_pos(probability_query_pos<=0.5);
                CC = rho.*log2(rho)-sum(probability_obj*ones(1,actual_inq).*rho_k.*log2(rho_k))';
                [minCC(nq), minCCid] = min(CC);    
                if minCCid(1)>inq/2
                    strategy(nq) = 1;
                else
                    strategy(nq) = 2;
                end
                
                query_id = aid(find(CC==min(CC),1));
                [options(1),options(2)] = find(queryID==query_id);

                if isempty(options)
                    probability_obj_set(:,nq:end) = repmat(probability_obj,1,MAX_ITER-nq+1);
                    expected_values(:,nq:end) = repmat(expected_value,1,MAX_ITER-nq+1);
                    partworths(XID,nq:end) = repmat(w0',1,MAX_ITER-nq+1);
                    nq = MAX_ITER+1;
                else
                    queryID(min(options),max(options))=0;
                    w_now = mvnrnd(w,Dw,1)'*theta;
                    u = (Xf(options(1),:)-Xf(options(2),:))*w_now;
                    if rand()<1/(1+exp(-u))
                        pairs(nq,:) = [options(1),options(2)];
                    else
                        pairs(nq,:) = [options(2),options(1)];
                    end

                    probability_obj_set(:,nq) = probability_obj;
                    partworths(XID,nq) = w0';
                    expected_values(:,nq) = expected_value;
                    nq = nq+1;
                    nquery = nquery - 1;
                end
                toc
                aaa=1;
            else
                probability_obj_set(:,nq:end) = repmat(probability_obj,1,MAX_ITER-nq+1);
                expected_values(:,nq:end) = repmat(expected_value,1,MAX_ITER-nq+1);
                partworths(XID,nq:end) = repmat(w0',1,MAX_ITER-nq+1);
                nq = MAX_ITER+1;
            end
        else
            options = randperm(nx,2); % random initial query
            w_now = mvnrnd(w,Dw,1)'*theta;
            u = (Xf(options(1),:)-Xf(options(2),:))*w_now;
            if rand()<1/(1+exp(-u))
                pairs(nq,:) = [options(1),options(2)];
            else
                pairs(nq,:) = [options(2),options(1)];
            end
            
            queryID(min(options),max(options))=0;
            probability_obj_set(:,nq) = probability_obj;
            expected_values(:,nq) = expected_value;
            partworths(XID,nq) = w0';
            nq = nq+1;
            nquery = nquery - 1;
        end
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    prob_set{test} = probability_obj_set;
    pairs_set{test} = pairs;
    partworths_set{test} = partworths;
    cond_set{test} = conds;
    expected_value_set{test} = expected_values;
    CC_set{test} = minCC;
    strategy_set{test} = strategy;
end
% save(['gisa_s',num2str(s),'_inq',num2str(inq),'_n',num2str(MAX_ITER),...
%     '_comp',num2str(num_competitor),'_theta',num2str(theta),...
%     '_nt',num2str(nt),'_09062016.mat'],...
%     'prob_set','pairs_set','partworths_set','target_best_set','cond_set',...
%     'expected_value_set','strategy_set','-v7.3');
% save(['profit_s',num2str(s),'_inq',num2str(inq),...
%     '_n',num2str(MAX_ITER),'_comp',num2str(num_competitor),...
%     '_theta',num2str(theta),'_0515.mat'],'prob_set','pairs_set',...
%     'partworths_set','target_best_set','strategy_set','-v7.3');