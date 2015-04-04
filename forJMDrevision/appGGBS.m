function [probability_obj_set,pairs] = ...
    appGGBS(s,d,Xf,dX,dXID,c,inq,const,alg,queryID,pairs,w,target,nq)
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
%             if nq>20
%                 if mean(probability_obj_set(1704,(nq-10):(nq-1)))<mean(probability_obj_set(1704,(nq-20):(nq-11)))
%                     stop = 1;
%                 end
%             end
%             non_zero_id = probability_obj>0;
%             non_zero = sum(non_zero_id);
            [~,sort_obj] = sort(probability_obj,'descend');
%             sort_obj = sort_obj(1:non_zero); %only need the non-zero ones

%             sort_obj = randperm(nx,inq);
            a = [];
            aid = [];
            count = 1;
%             outcount = 1;
            
            for i = 1:(length(sort_obj)-1)
                if count > inq
                    break;
                end
                for j = (i+1):length(sort_obj)
                    if count > inq
                        break;
                    end
%                     outcount = outcount + 1;
                    ii = min([sort_obj(i),sort_obj(j)]);
                    II = max([sort_obj(i),sort_obj(j)]);
%                     if queryID(ii,II)>0 && probability_obj(ii)>0 && probability_obj(II)>0
                    if queryID(ii,II)>0
                        a = [a;[ii,II]];
                        aid = [aid;queryID(ii,II)];
                        count = count + 1;
                    end
                end
            end
            
            actual_inq = length(aid);
            if actual_inq>0
                temp_set = zeros(nx,actual_inq);

    %             Q_remain = dX(queryList>0,:); % all remaining queries
                Q_remain = dX(aid,:);
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
    %             query_id = find(queryList>0,find(C==min(C),1));
    %             [options(1),options(2)] = find(queryID==query_id);
                query_id = aid(find(C==min(C),1));
                options = a(find(C==min(C),1),:);

                if isempty(options)
                    probability_obj_set = [probability_obj_set, probability_obj];
                    fprintf('no feasible query. terminated.');
                    return;
                end
    %             options = options(1,:);
                queryID(min(options),max(options))=0;
                queryList(queryID==query_id)=0;

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