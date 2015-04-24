clear; close all;
rng(6);
% 2D example
d = 2;
nx = 1000;

% Xf = [-1,-1;-1,1;1,-1;1,1];
Xf = rand(nx,d)*2-1;
% Xf = [cos(0:(pi/20):(2*pi));sin(0:(pi/20):(2*pi))]'+2;
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
% c = [1;1;10;10]/100;
% c = Xf(:,2);
c = ceil(rand(nx,1)*25)/100; %price - cost
% w = [1,-0.6];
w = rand(1,d)*2-1;

% available questions 
% Q = [cos(0:(pi/20):(2*pi));sin(0:(pi/20):(2*pi))]'+2;
% nq = size(Q,1); % number of queries
addpath('..\..\Tools\liblinear\matlab');


TEST = 1;
MAX_ITER = 200;
prob_set = cell(TEST,1);
pairs_set = cell(TEST,1);
partworths_set = cell(TEST,1);
s = 1e3;
inq = min(10,nq);
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

fprintf('truth: %d',target_best);

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
    guesses = zeros(MAX_ITER,1);
    nq = 1;
    nquery = sum(sum(queryID>0));
    queryList = ones(nquery,1);
%     A = zeros(MAX_ITER,d); % question set
    
    % to plot everything
    figure; hold on;
    xrange = -10:1:30;
    yrange = -10:1:30;
    [X,Y] = meshgrid(xrange,yrange);
    Z = zeros(length(yrange),length(xrange));
    ZZ = zeros(length(yrange),length(xrange));
    for i = 1:length(xrange)
        for j = 1:length(yrange)
            xx = [xrange(i);yrange(j)];
            util = Xf(1:nt,:)*xx;
            obj_app = bsxfun(@plus,-log(1+exp(-util)),log(c(1:nt,:)));
            [~,best_temp] = max(obj_app);
            Z(j,i) = best_temp;
            if best_temp==target_best
                Z(j,i) = nx*2;
            end
        end
    end
    
    % find the next query
    while  nq>=1
      
        % calculate A
        A = dX(DX(sub2ind(size(DX),pairs(1:nq-1,1),pairs(1:nq-1,2))),:);
        A = bsxfun(@times,A,sign(pairs(1:nq-1,2)-pairs(1:nq-1,1)));
        [probability_obj,W,I,unique_I,w0,C] = appObjDistribution(s,d,W0,Xf(1:nt,:),c(1:nt,:),A(1:nq-1,:));
        
%         util = Xf(1:nt,:)*w0';
%         obj_app = bsxfun(@plus,-log(1+exp(-util)),log(c(1:nt,:)));
%         [~,guess] = max(obj_app);
        [~,guess] = max(probability_obj);
        fprintf('iter: %d, truth: %d, guess: %d \n',nq, target_best, guess);
        
        if nq>MAX_ITER
            figure;
            plot(guesses);   
%             probability_obj_set = [probability_obj_set, probability_obj];
%             partworths = [partworths, w0'];
%             fprintf('max query number exceeded. terminated.');
            return;
        else
            [~,sort_obj] = sort(probability_obj,'descend');
            sort_obj_nonzero = sort_obj(1:sum(probability_obj>0));
            a = [];
            aid = [];
            count = 1;
            
            for i = 1:(length(sort_obj_nonzero)-1)
                if count > inq
                    break;
                end
                for j = (i+1):length(sort_obj_nonzero)
                    if count > inq
                        break;
                    end
                    ii = min([sort_obj_nonzero(i),sort_obj_nonzero(j)]);
                    II = max([sort_obj_nonzero(i),sort_obj_nonzero(j)]);
                    if queryID(ii,II)>0
                        a = [a;[ii,II]];
                        aid = [aid;queryID(ii,II)];
                        count = count + 1;
                    end
                end
            end
            
            if count <= inq
                sort_obj_zero = setdiff((1:nx)',sort_obj_nonzero);
                sort_obj_zero = sort_obj_zero(randperm(length(sort_obj_zero)));
                for i = 1:length(sort_obj_nonzero)
                    if count > inq
                        break;
                    end
                    for j = 1:length(sort_obj_zero)
                        if count > inq
                            break;
                        end
                        ii = min([sort_obj_nonzero(i),sort_obj_zero(j)]);
                        II = max([sort_obj_nonzero(i),sort_obj_zero(j)]);
                        if queryID(ii,II)>0
                            a = [a;[ii,II]];
                            aid = [aid;queryID(ii,II)];
                            count = count + 1;
                        end
                    end
                end
            end
            
%             actual_inq = inq;
            actual_inq = length(aid);
            if actual_inq>0
                query_id = aid;
                options = a;
                temp_set = zeros(nt,actual_inq);
                Q_remain = dX(aid,:);
%                 Q_remain = Q;
                u = W*Q_remain'>0; % get all utility signs
                u = bsxfun(@times, u, I); % label all positive signs with query numbers
                for j = 1:length(unique_I)
                    temp_set(unique_I(j),:) = sum(u==unique_I(j),1)/sum(I==unique_I(j));
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
                
%                 if sum(CC.^2)<1e-16 % all CC = 0
%                     A = dX(DX(sub2ind(size(DX),pairs(1:nq-1,1),pairs(1:nq-1,2))),:);
%                     A = bsxfun(@times,A,sign(pairs(1:nq-1,2)-pairs(1:nq-1,1)));
% %                     model = train(ones(size(A,1),1),...
% %                         sparse(A), sprintf('-s 0 -e 1e-6 -q -c %f', C));
% %                     w0 = (model.w)';
% %                     [w0, C] = crossvalidation(sparse(A),ones(size(A,1),1));
% %                     w0 = (A'*A+1/C*eye(size(A,2)))\A'*ones(size(A,1),1);
%                     B = (eye(size(A,2))-(w0'*w0)/(w0*w0'))*(A'*A+1/C*eye(size(A,2)));
%                     [V,D] = eig(B);
%                     e = diag(D);
%                     v = V(:,e==min(e(e>1e-3)));
%                     v = v(:,1);
%                     s1 = dX*w0';
%                     s2 = (dX*v)./sqrt(sum(dX.^2,2));
%                     unsampled = unique(queryID);
%                     unsampled(unsampled==0)=[];
%                     s1_id = find(abs(s1(unsampled))==min(abs(s1(unsampled))));
%                     s2_id = find(abs(s2(unsampled(s1_id)))==max(abs(s2(unsampled(s1_id)))));
%                     query_id = unsampled(s1_id(s2_id(1)));
%                     [a1,a2] = find(dXID==query_id);
%                     options = [a1,a2];
%                 else
                    query_id = aid(find(CC==min(CC),1));
                    options = a(find(CC==min(CC),1),:);
%                 end
%                 options = Q(find(CC==min(CC),1),:);
                if isempty(options)
                    probability_obj_set(:,nq) = probability_obj;
                    partworths(:,nq) = mean(w0,1)';
                    fprintf('no feasible query. terminated.');
                    return;
                end
%                 queryID(min(options),max(options))=0;
%                 queryList(queryID==query_id)=0;
                
                w_now = mvnrnd(w,Dw,1)'*theta;
                u = (Xf(options(1),:)-Xf(options(2),:))*w_now;
%                 u = options*w_now;
                if rand()<1/(1+exp(-u))
                    pairs(nq,:) = [options(1),options(2)];
%                     A(nq,:) = options;
                else
                    pairs(nq,:) = [options(2),options(1)];
%                     A(nq,:) = -options;
                end

                probability_obj_set(:,nq) = probability_obj;
                partworths(:,nq) = mean(w0,1)';
                guesses(nq) = guess==target_best;
                nq = nq+1;
                nquery = nquery - 1;
                
                
                %%%%%%%%%%%%%%%%%%%%%%%% plot things %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                cla; hold on;
                
%                 A = dX(DX(sub2ind(size(DX),pairs(1:(nq-1),1),pairs(1:(nq-1),2))),:);
%                 A = bsxfun(@times,A,sign(pairs(1:(nq-1),2)-pairs(1:(nq-1),1)));
                contourf(X,Y,Z);colormap summer;
                

                
                for Cid = 1:length(C)
%                     As = bsxfun(@times, A, sqrt(exp(A*w0(Cid,:)'))./(1+exp(A*w0(Cid,:)')));
%                     As(isnan(As))=0;
%                     Sigma_inv = (eye(d)/C(Cid)+(As'*As));
% 
%                     for i = 1:length(xrange)
%                         for j = 1:length(yrange)
%                             xx = [xrange(i);yrange(j)];
%                             ZZ(j,i) = exp(-1/2*(xx'-w0(Cid,:))*Sigma_inv*(xx-w0(Cid,:)'));
%                         end
%                     end
%                     contour(X,Y,ZZ); 
                    plot(w0(Cid,1),w0(Cid,2),'.b','MarkerSize',30)
                end
                
                
                plot(wtrue(1),wtrue(2),'.r','MarkerSize',50)
                plot(Xf(pairs(nq-1,1),1),Xf(pairs(nq-1,1),2),'sb','MarkerSize',15);
                plot(Xf(pairs(nq-1,2),1),Xf(pairs(nq-1,2),2),'sb','MarkerSize',15);
                plot(Xf(target_best,1),Xf(target_best,2),'sy','MarkerSize',15,'MarkerFaceColor','y');
                plot(Xf(guess,1),Xf(guess,2),'sr','MarkerSize',10,'MarkerFaceColor','r');
                for i = 1:nx
                    plot(Xf(i,1),Xf(i,2),'sk','MarkerSize',max(probability_obj(i)*80,5));
                end
                
                if nq>2
                    k = -(Xf(pairs(nq-2,1),1)-Xf(pairs(nq-2,2),1))/(Xf(pairs(nq-2,1),2)-Xf(pairs(nq-2,2),2));                    k = -(Xf(pairs(nq-2,1),1)-Xf(pairs(nq-2,2),1))/(Xf(pairs(nq-2,1),2)-Xf(pairs(nq-2,2),2));
%                     k = -A(nq-2,1)/A(nq-2,2);
                    plot(xrange,k*xrange,'b','linewidth',2);
                end
                %%%%%%%%%%%%%%%%%%%%%%%% plot things %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                wait = 1;
                
            else
                probability_obj_set(:,nq) = probability_obj;
                partworths(:,nq) = mean(w0,1)';
                guesses(nq) = guess==target_best;
                fprintf('no feasible query. terminated.');
                return;
            end
        end
    end    
    
    
    
    
    prob_set{test} = probability_obj_set;
    pairs_set{test} = pairs;
    partworths_set{test} = partworths;
end

