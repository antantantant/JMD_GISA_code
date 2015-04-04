function [probability_obj_set,pairs] = ...
    halfV(s,d,Xf,dX,dXID,c,inq,const,alg,queryID,pairs,w,target,nq)
    nx = size(Xf,1);
    
    % calculate A
    alg = 'profit';
    [A,W] = appObjDistribution(s,d,w,Xf,dX,dXID,c,const,alg,pairs,[]);
    probability_obj = A/sum(A);
    fprintf('|%d, %f|',nq,JSdivergence(probability_obj,target));
    if isnan(JSdivergence(probability_obj,target))
        what = 1;
    end
    alg = 'pref_only';
    
    % find the next query
    if  nq>=0
%         if ~isempty(pairs)&&sum(probability_obj<(1/nx*1e-3))==nx-1
        if nq>100
            probability_obj_set = probability_obj;
            return;
        else
            if ~isempty(pairs)
                DX = dXID + dXID';
                A = dX(DX(sub2ind(size(DX),pairs(:,1),pairs(:,2))),:);
                A = bsxfun(@times,A,sign(pairs(:,2)-pairs(:,1)));
%                 model = train([ones(size(A,1),1);-ones(size(A,1),1)],...
%                         sparse([A;-A]),'-s 0 -c 1e6 -e 1e-6 -q');
%                 w0 = model.w;
%                 corr = dX*w0';
                
%                 zerorow = find(sum(dX.^2,2)<1e-6); %find all-zero rows
                w0 = (A'*A+1/(nq-1)*eye(size(A,2)))\A'*ones(size(A,1),1);
                B = (eye(size(A,2))-w0*w0'/(w0'*w0))*(A'*A+1/(nq-1)*eye(size(A,2)));
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
%                 [~,sort_obj] = sort(probability_obj,'descend');
%                 a = sort_obj(1:min(inq,nx));
%                 aa = queryID(a,a);
%                 aaa = unique(reshape(aa,1,size(aa,1)*size(aa,2)));
%                 aaa(aaa==0)=[];
%                 a = zeros(length(aaa),2);
%                 for i = 1:length(aaa)
%                     [a(i,1),a(i,2)] = find(queryID==aaa(i));
%                 end
%                 no = nx;
%                 temp_set = zeros(no,size(a,1));
%                 for i = 1:size(a,1)
%                     % calculate the probability of a new query returning 1
%                     temp_set(:,i) = appObjDistribution(s,d,W,Xf,dX,dXID,c,const,alg,[pairs;a(i,:)],[]);
%                 end
%                 probability_query_pos = sum(temp_set,1)';% probability of new query being positive
%                 cut = isnan(probability_query_pos);
%                 if sum(cut)>0
%                     a(cut,:) = [];
%                     probability_query_pos(cut)=[];
%                     temp_set(:,cut)=[];
%                 end
%                 rho = zeros(size(a,1),1);
%                 probability_query_pos(probability_query_pos==0)=1e-99;
% 
%                 rho(probability_query_pos>0.5) = probability_query_pos(probability_query_pos>0.5);
%                 rho(probability_query_pos<=0.5) = 1-probability_query_pos(probability_query_pos<=0.5);
%                 C = rho.*log2(rho);
%                 options = a(C==min(C),:);
                options = [2360,85];
            end
            
            if isempty(options)
                probability_obj_set = probability_obj;
                return;
            end
            options = options(1,:);
            queryID(min(options),max(options))=0;
            
            u = (Xf(options(1),:)-Xf(options(2),:))*w;
            if u>0
                pairs = [pairs;options(1),options(2)];
            elseif u<0
                pairs = [pairs;options(2),options(1)];
            end
            [probability_obj_sub_set,pairs] = halfV(s,d,Xf,dX,dXID,c,inq,const,alg,queryID,pairs,w,target,nq+1);
            probability_obj_set = [probability_obj, probability_obj_sub_set];
        end
    end