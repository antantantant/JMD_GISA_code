function [probability_obj_set,pairs] = ...
    appGGBS(s,d,Xf,dX,dXID,c,inq,const,alg,queryID,pairs,w,target,nq)
    nx = size(Xf,1);
    
    % calculate A
    [A,W] = appObjDistribution(s,d,w,Xf,dX,dXID,c,const,alg,pairs,[]);
    probability_obj = A/sum(A);
    fprintf('|%d, %.2f|',nq,JSdivergence(probability_obj,target));
    % find the next query
    if  nq>=0
%         if ~isempty(pairs)&&sum(probability_obj<(1/nx*1e-3))==nx-1
        if nq>100
            probability_obj_set = probability_obj;
            return;
        else
            [~,sort_obj] = sort(probability_obj,'descend');
            a = sort_obj(1:min(inq,nx));
            aa = queryID(a,a);
            aaa = unique(reshape(aa,1,size(aa,1)*size(aa,2)));
            aaa(aaa==0)=[];
            a = zeros(length(aaa),2);
            for i = 1:length(aaa)
                [a(i,1),a(i,2)] = find(queryID==aaa(i));
            end
            no = nx;
            temp_set = zeros(no,size(a,1));
            for i = 1:size(a,1)
                % calculate the probability of a new query returning 1
                temp_set(:,i) = appObjDistribution(s,d,W,Xf,dX,dXID,c,const,alg,[pairs;a(i,:)],[]);
%                 if sum(temp_set(:,i)>0)==0
%                      probability_obj_set = [probability_obj,zeros(size(probability_obj))];
%                      return;
%                 end
            end
            probability_query_pos = sum(temp_set,1)';% probability of new query being positive
            cut = isnan(probability_query_pos);
            if sum(cut)>0
                a(cut,:) = [];
                probability_query_pos(cut)=[];
                temp_set(:,cut)=[];
            end
            rho_k = zeros(no,size(a,1));
            rho = zeros(size(a,1),1);
            probability_query_pos(probability_query_pos==0)=1e-99;

            temp = bsxfun(@rdivide,temp_set,probability_obj+1e-99);
            rho_k(temp>0.5) = temp(temp>0.5);
            rho_k(temp<=0.5) = 1-temp(temp<=0.5);
            rho(probability_query_pos>0.5) = probability_query_pos(probability_query_pos>0.5);
            rho(probability_query_pos<=0.5) = 1-probability_query_pos(probability_query_pos<=0.5);
            C = rho.*log2(rho)-sum(probability_obj*ones(1,size(a,1)).*rho_k.*log2(rho_k))';
            options = a(C==min(C),:);

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
            [probability_obj_sub_set,pairs] = appGGBS(s,d,Xf,dX,dXID,c,inq,const,alg,queryID,pairs,w,target,nq+1);
            probability_obj_set = [probability_obj, probability_obj_sub_set];
        end
    end