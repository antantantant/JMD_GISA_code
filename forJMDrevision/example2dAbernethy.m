clear;
rng(0);
% 2D example
nx = 1000;
d = 3;
Xf = rand(nx,d)*2-1;
nq = nx*(nx-1)/2; % number of queries
dX = zeros(nq,d); % attribute difference data
dXID = zeros(nx,nx); % pair to query id mapping
count = 1;
for i = 1:nx
    for j = i+1:nx
        dXID(i,j) = count;
        dX(count,:) = Xf(i,:) - Xf(j,:);
        count = count + 1;
    end
end
c = rand(nx,1); %price - cost
w = rand(1,d)*2-1;

addpath('..\..\Tools\liblinear\matlab');


TEST = 1;
MAX_ITER = min(200,nq);
prob_set = cell(TEST,1);
pairs_set = cell(TEST,1);
partworths_set = cell(TEST,1);
s = 1e3;
inq = min(nq,10);
W0 = mvnrnd(zeros(1,d),eye(length(w)),s);
alg = '';
sigma = 1e-6;
Dw = eye(d)*sigma; % randomness in user choices

nt = size(Xf,1); % number of testing object

% Calculate true distribution under sigma
Wtrue = mvnrnd(w,Dw,s);
theta = 10;
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
    partworths = zeros(d,MAX_ITER);
    pairs = zeros(MAX_ITER,2);
    nq = 1;
    nquery = sum(sum(queryID>0));
    queryList = ones(nquery,1);
        
    % find the next query
    while  nq>=1
        C = nq; % update C
        
        % calculate A
        [probability_obj,~,~,~,w0,C] = appObjDistribution(s,d,W0,Xf(1:nt,:),dX,dXID,c(1:nt,:),alg,pairs(1:nq-1,:));
        w0 = w0';
        fprintf('%d, %.2f, %.2f %.2f \n',nq,...
            probability_obj(target_best), corr(w0,wtrue), norm(w0-wtrue));

        if nq>MAX_ITER
            return;
        else
            if nq>1
                
                A = dX(DX(sub2ind(size(DX),pairs(1:nq-1,1),pairs(1:nq-1,2))),:);
                A = bsxfun(@times,A,sign(pairs(1:nq-1,2)-pairs(1:nq-1,1)));
%                 model = train(ones(size(A,1),1),...
%                             sparse(A), sprintf('-s 0 -e 1e-6 -q -c %f', C));
%                 w0 = (model.w)';
%                 [w0, C] = crossvalidation(sparse(A),ones(size(A,1),1));
               
%                 w0 = (A'*A+1/C*eye(size(A,2)))\A'*ones(size(A,1),1);
                B = (eye(size(A,2))-w0*w0'/(w0'*w0))*(A'*A+1/C*eye(size(A,2)));
                [V,D] = eig((B+B')/2);
                e = diag(D);
                v = V(:,e==min(e(e>1e-12)));
                v = v(:,1);
                s1 = dX*w0;
                s2 = (dX*v)./sqrt(sum(dX.^2,2));
                                
                unsampled = unique(queryID);
                unsampled(unsampled==0)=[];
%                 unsampled = setdiff(unsampled, zerorow);
                s1_id = find(abs(s1(unsampled))==min(abs(s1(unsampled))));
                s2_id = find(abs(s2(unsampled(s1_id)))==max(abs(s2(unsampled(s1_id)))));
                id = unsampled(s1_id(s2_id(1)));
                [a1,a2] = find(dXID==id);
                options = [a1,a2];
            else
                options = randperm(nx,2); % random initial query
            end
            
            if isempty(options)
                probability_obj_set(:,nq) = probability_obj;
                partworths(:,nq) = w0;
                return;
            end
            options = options(1,:);
            queryID(min(options),max(options))=0;
            
            
            w_now = mvnrnd(w,Dw,1)'*theta;
            u = (Xf(options(1),:)-Xf(options(2),:))*w_now;
            if rand()<1/(1+exp(-u))
                pairs(nq,:) = [options(1),options(2)];
            else
                pairs(nq,:) = [options(2),options(1)];
            end
            
            probability_obj_set(:,nq) = probability_obj;
            partworths(:,nq) = w0;
            nq = nq+1;
        end
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    plot(probability_obj_set(target_best,:));
    
    
    
    prob_set{test} = probability_obj_set;
    pairs_set{test} = pairs;
    partworths_set{test} = partworths;
end