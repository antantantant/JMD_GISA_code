function [probability_obj_set,pairs,partworths] = ...
    halfV(s,d,Xf,dX,dXID,c,alg,queryID,W0,w,Dw,target,C)
    nx = size(Xf,1);
    DX = dXID + dXID';
    probability_obj_set = [];
    partworths = [];
    pairs = [];
    nq = 0;
    w0 = zeros(d,1);
    % find the next query
    while  nq>=0
        % calculate A
%         A = appObjDistribution(s,d,W0,Xf,dX,dXID,c,alg,pairs,C);
%         probability_obj = A/sum(A);
%         fprintf('%d, %.2f, %.2f %.2f \n',nq,JSdivergence(probability_obj,target),...
%             probability_obj(1704), corr(w0,w));

%         if ~isempty(pairs)&&sum(probability_obj<(1/nx*1e-3))==nx-1
        if nq>500
%             probability_obj_set = [probability_obj_set, probability_obj];
            partworths = [partworths, w0];
            return;
        else
            if ~isempty(pairs)
                C = (nq+1); % update C
                A = dX(DX(sub2ind(size(DX),pairs(:,1),pairs(:,2))),:);
                A = bsxfun(@times,A,sign(pairs(:,2)-pairs(:,1)));
%                 model = train([ones(size(A,1),1);-ones(size(A,1),1)],...
%                         sparse([A;-A]),'-s 0 -c 1e6 -e 1e-6 -q');
%                 w0 = model.w;
%                 corr = dX*w0';
                
                w0 = (A'*A+1/C*eye(size(A,2)))\A'*ones(size(A,1),1);
                B = (eye(size(A,2))-w0*w0'/(w0'*w0))*(A'*A+1/C*eye(size(A,2)));
                [V,D] = eig(B);
                e = diag(D);
                v = V(:,find(e==min(e(e>1e-3))));
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
%                 probability_obj_set = [probability_obj_set, probability_obj];
                partworths = [partworths, w0];
                return;
            end
            options = options(1,:);
            queryID(min(options),max(options))=0;
            
            
            w_true = mvnrnd(w,Dw,1);
            u = (Xf(options(1),:)-Xf(options(2),:))*w_true';
            if u>0
                pairs = [pairs;options(1),options(2)];
            elseif u<0
                pairs = [pairs;options(2),options(1)];
            end
            
%             probability_obj_set = [probability_obj_set, probability_obj];
            partworths = [partworths, w0];
            nq = nq+1;
        end
    end