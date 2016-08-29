function [probability_obj_set,pairs,partworths] = ...
    halfV(s,d,Xf,dX,dXID,c,alg,queryID,W0,w,Dw,target,target_best,MAX_ITER)
    nx = size(Xf,1);
    DX = dXID + dXID';
    probability_obj_set = zeros(10,MAX_ITER);
    partworths = zeros(d,MAX_ITER);
    pairs = zeros(MAX_ITER,2);
    nq = 1;
    w0 = zeros(d,1);
    
    % find the next query
    while  nq>=1
%         C = nq; % update C
        
        % calculate A
        A = appObjDistribution(s,d,W0,Xf(1:10,:),dX,dXID,c(1:10,:),alg,pairs(1:nq-1,:),C);
        probability_obj = A/sum(A);
        fprintf('%d, %.2f, %.2f %.2f \n',nq,JSdivergence(probability_obj,target),...
            probability_obj(target_best), corr(w0,w));

        if nq>MAX_ITER
            return;
        else
            if nq>1
                
                A = dX(DX(sub2ind(size(DX),pairs(1:nq-1,1),pairs(1:nq-1,2))),:);
                A = bsxfun(@times,A,sign(pairs(1:nq-1,2)-pairs(1:nq-1,1)));
%                 model = train(ones(size(A,1),1),...
%                             sparse(A), sprintf('-s 0 -e 1e-6 -q -c %f', C));
%                 w0 = (model.w)';
                [w0, C] = crossvalidation(sparse(A),ones(size(A,1),1));
                
%                 % normalize w according to the original scale paper
%                 w0(1:5) = w0(1:5)/abs(sum(w0(1:5)))*0.02;
%                 w0(6:10) = w0(6:10)/abs(sum(w0(6:10)))*0.02;
%                 w0(11:15) = w0(11:15)/abs(sum(w0(11:15)))*0.02;
%                 w0(16:20) = w0(16:20)/abs(sum(w0(16:20)))*0.02;
%                 w0(21:25) = w0(21:25)/abs(sum(w0(21:25)))*0.02;
%                 w0(26:30) = w0(26:30)/abs(sum(w0(26:30)))*0.02;
                
%                 w0 = (A'*A+1/C*eye(size(A,2)))\A'*ones(size(A,1),1);
                B = (eye(size(A,2))-w0*w0'/(w0'*w0))*(A'*A+1/C*eye(size(A,2)));
                [V,D] = eig(B);
                e = diag(D);
                v = V(:,e==min(e(e>1e-3)));
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
            
            
            w_true = mvnrnd(w,Dw,1);
            u = (Xf(options(1),:)-Xf(options(2),:))*w_true';
            if u>0
                pairs(nq,:) = [options(1),options(2)];
            elseif u<0
                pairs(nq,:) = [options(2),options(1)];
            end
            
            probability_obj_set(:,nq) = probability_obj;
            partworths(:,nq) = w0;
            nq = nq+1;
        end
    end