function [probability_obj_set,pairs,method] = ...
    appGGBS_PLUS(s,d,Xf,dX,dXID,c,inq,const,alg,queryID,pairs,w,target,nq)
    nx = size(Xf,1);
    DX = dXID + dXID';
    probability_obj_set = [];
    nquery = sum(sum(queryID>0));
    queryList = ones(nquery,1);
    method = zeros(10001,1);
    
    % find the next query
    while  nq>=0
        % calculate A
        [A,W,I,unique_I] = appObjDistribution(s,d,w,Xf,dX,dXID,c,const,alg,pairs,[]);
        probability_obj = A/sum(A);
        fprintf('%d, %.2f, %.2f \n',nq,JSdivergence(probability_obj,target),...
            probability_obj(1704));
        
%         if ~isempty(pairs)&&sum(probability_obj<(1/nx*1e-3))==nx-1
        if nq>1000
            probability_obj_set = [probability_obj_set, probability_obj];
            fprintf('max query number exceeded. terminated.');
            return;
        else
            if nq>200
                if mean(probability_obj_set(1704,(nq-100):(nq-1)))<mean(probability_obj_set(1704,(nq-200):(nq-101)))
                    stop = 1;
                end
            end
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
            
            if count < inq
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
            
            
            actual_inq = length(aid);
            if actual_inq>0
                temp_set = zeros(nx,actual_inq);
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
                
                if sum(C.^2)<1e-16 % all C = 0
                    A = dX(DX(sub2ind(size(DX),pairs(:,1),pairs(:,2))),:);
                    A = bsxfun(@times,A,sign(pairs(:,2)-pairs(:,1)));
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
                    s1_id = find(abs(s1(unsampled))==min(abs(s1(unsampled))));
                    s2_id = find(abs(s2(unsampled(s1_id)))==max(abs(s2(unsampled(s1_id)))));
                    id = unsampled(s1_id(s2_id(1)));
                    [a1,a2] = find(dXID==id);
                    options = [a1,a2];
                    
                    method(nq) = 1;
                else
                    query_id = aid(find(C==min(C),1));
                    options = a(find(C==min(C),1),:);
                    
                    method(nq) = 0;
                end

                if isempty(options)
                    probability_obj_set = [probability_obj_set, probability_obj];
                    fprintf('no feasible query. terminated.');
                    return;
                end
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