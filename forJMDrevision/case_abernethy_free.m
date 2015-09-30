% matlabpool open;
% load('..\basedata.mat');
addpath('..\..\\Tools\liblinear\matlab');

c = Xf(:,26:30)*price'-cv; %price - cost
TEST = 12;
MAX_ITER = 20000;
prob_set = cell(TEST,1);
dxs_set = cell(TEST,1);
partworths_set = cell(TEST,1);
target_best_set = zeros(TEST,1);
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

theta = 1;
wtrue = w*theta;
wtrue(1:5) = wtrue(1:5)-wtrue(5);
wtrue(6:10) = wtrue(6:10)-wtrue(10);
wtrue(11:15) = wtrue(11:15)-wtrue(15);
wtrue(16:20) = wtrue(16:20)-wtrue(20);
wtrue(21:25) = wtrue(21:25)-wtrue(25);
wtrue(26:30) = wtrue(26:30)-wtrue(30);

num_competitor = 1;

parfor test = 1:TEST
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
    exp_sum_delta_util = exp_delta_util*repmat(eye(nt),num_competitor,1)+1;
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
    dxs = zeros(MAX_ITER,length(XID));
    conds = zeros(1,MAX_ITER);
%     probability_obj_set = [prob_set{test},zeros(2455,500)];
%     partworths = [partworths_set{test},zeros(30,500)];
%     pairs = [pairs_set{test};zeros(500,2)];
    nq = 1;
    nquery = sum(sum(queryID>0));
    queryList = ones(nquery,1);
        
    % find the next query
    while  nq<=MAX_ITER
        % calculate A
        A = dxs(1:nq-1,:);
        [probability_obj,W,I,unique_I,w0,C,weights] = ...
            appObjDistribution(s,length(XID),W0(:,XID),Xf(1:nt,XID),c(1:nt,:),A,competitors);
        util = Xf(1:nt,XID)*w0';
        obj_app = bsxfun(@plus,-log(1+exp(-util)),log(c(1:nt,:)));
%         
%         
%         % test corr only
%         [w0, C, weights] = crossvalidation(sparse(A),ones(size(A,1),1)); 
%         probability_obj = zeros(nt,1);
%         if isempty(w0) 
%             w0 = zeros(1,length(XID));
%         end
%         % test corr only
            
        if nq>1
            As = bsxfun(@times, A, sqrt(exp(A*w0'))./(1+exp(A*w0')));
            As(isnan(As))=0;
            Sigma_inv = (eye(length(XID))/C+(As'*As));
            conds(nq) = cond(Sigma_inv);
            
            [~,guess] = max(probability_obj);
            fprintf('iter: %d, truth: %d (%f), guess: %d, point guess: %d, corr: %0.2f, norm: %0.2f  \n',...
                nq, target_best, probability_obj(target_best), guess,...
                find(obj_app==max(obj_app),1), corr(w0',wtrue(XID)'), norm(w0));
            if corr(w0',wtrue(XID)')>0.9
                wait = 1;
            end
            
%             fprintf('iter: %d, corr: %0.2f, cond: %0.2f \n\n\n\n',...
%                 nq, corr(w0',wtrue(XID)'), conds(nq));


            % Abernethy's method, do not project onto the subspace
            % perpendicular to w0
            B = eye(size(A,2))/C+(As'*As);
            % Abernethy's method, w/ logistic hessian
%             B = (eye(size(A,2))-(w0'*w0)/(w0*w0'))*(eye(size(A,2))/C+(As'*As));
            % Abernethy's method, w/ quadratic hessian
%             B = (eye(size(A,2))-(w0'*w0)/(w0*w0'))*(A'*A+1/C*eye(size(A,2)));
            [V,D] = eig(B);
            e = diag(D);
            v =  V(:,e==min(e(e>1e-12)));
            v = v(:,1);
            
%             options = v'/norm(v)-w0/norm(w0);
%             options = options/norm(options);
            options = v'/norm(v);
        else
            options = rand(1,length(XID)); % random initial dx
        end

        if isempty(options)
            probability_obj_set(:,nq:end) = repmat(probability_obj,1,MAX_ITER-nq+1);
            partworths(XID,nq:end) = repmat(w0',1,MAX_ITER-nq+1);
            nq = MAX_ITER+1;
        else
%             w_now = mvnrnd(wtrue,Dw,1)'*theta;
            w_now = wtrue(XID)'*theta;
            u = options*w_now;
            if rand()<1/(1+exp(-u))
                dxs(nq,:) = options;
            else
                dxs(nq,:) = -options;
            end

            probability_obj_set(:,nq) = probability_obj;
            partworths(XID,nq) = w0';
            nq = nq+1;               
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    prob_set{test} = probability_obj_set;
    pairs_set{test} = pairs;
    partworths_set{test} = partworths;
    cond_set{test} = conds;
    dxs_set{test} = dxs;
end
save(['abernethy_free_s',num2str(s),'_n',num2str(MAX_ITER),...
    '_comp',num2str(num_competitor),'_theta',num2str(theta),'_0923.mat'],...
    'prob_set','dxs_set','partworths_set','target_best_set','cond_set','-v7.3');
