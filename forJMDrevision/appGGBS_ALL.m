function [probability_obj_set,pairs] = ...
    appGGBS_ALL(s,d,Xf,dX,dXID,c,inq,const,alg,queryID,pairs,w,target,nq)
    nx = size(Xf,1);
    probability_obj_set = [];
    nquery = sum(sum(queryID>0));
    queryList = ones(nquery,1);
    
    % find the next query
    while  nq>=0
        % calculate A
        [A,W,I,unique_I] = appObjDistribution(s,d,w,Xf,dX,dXID,c,const,alg,pairs,[]);
        probability_obj = A/sum(A);
        fprintf('%d, %.2f, %.2f \n',nq,JSdivergence(probability_obj,target),...
            probability_obj(1704));
        
%         if ~isempty(pairs)&&sum(probability_obj<(1/nx*1e-3))==nx-1
        if nq>10000
            probability_obj_set = [probability_obj_set, probability_obj];
            fprintf('max query number exceeded. terminated.');
            return;
        else
            actual_inq = sum(queryList);
            if actual_inq>0
                temp_set = zeros(nx,actual_inq);

                Q_remain = dX(queryList>0,:); % all remaining queries
                u = W*Q_remain'>0; % get all utility signs
                u = bsxfun(@times, u, I); % label all positive signs with query numbers
                for j = 1:length(unique_I)
                    temp_set(unique_I(j),:) = sum(u==unique_I(j),1)/sum(I==unique_I(j));
                end

                probability_query_pos = sum(bsxfun(@times,temp_set,probability_obj),1)';% probability of new query being positive

                rho_k = zeros(nx,actual_inq);
                rho = zeros(actual_inq,1);
                probability_query_pos(probability_query_pos==0)=1e-99;

                rho_k(temp_set>0.5) = temp_set(temp_set>0.5);
                rho_k(temp_set<=0.5) = 1-temp_set(temp_set<=0.5);
                rho(probability_query_pos>0.5) = probability_query_pos(probability_query_pos>0.5);
                rho(probability_query_pos<=0.5) = 1-probability_query_pos(probability_query_pos<=0.5);
                C = rho.*log2(rho)-sum(probability_obj*ones(1,actual_inq).*rho_k.*log2(rho_k))';
                remaining_query_id = find(queryList>0);
                query_id = remaining_query_id(find(C==min(C),1));
                [options(1), options(2)] = find(queryID==query_id);

                if isempty(options)
                    probability_obj_set = [probability_obj_set, probability_obj];
                    fprintf('no feasible query. terminated.');
                    return;
                end
    %             options = options(1,:);
                queryID(min(options),max(options))=0;
                queryList(query_id)=0;

                u = (Xf(options(1),:)-Xf(options(2),:))*w;
                if u>0
                    pairs = [pairs;options(1),options(2)];
                elseif u<0
                    pairs = [pairs;options(2),options(1)];
                end
                probability_obj_set = [probability_obj_set, probability_obj];
                nq = nq+1;
                nquery = nquery - 1;
            else
                probability_obj_set = [probability_obj_set, probability_obj];
                fprintf('no feasible query. terminated.');
                return;
            end
        end
    end