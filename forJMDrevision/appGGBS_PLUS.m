function [probability_obj_set,pairs,partworths] = ...
    appGGBS_PLUS(s,d,Xf,dX,dXID,c,inq,alg,queryID,W0,w,Dw,target,target_best,...
    MAX_ITER,...
    prob_set, pairs_set, partworths_set)
    nx = size(Xf,1);
    DX = dXID + dXID';
    probability_obj_set = zeros(nx,MAX_ITER);
    partworths = zeros(d,MAX_ITER);
    pairs = zeros(MAX_ITER,2);
    nq = 1;
    nquery = sum(sum(queryID>0));
    queryList = ones(nquery,1);
    
%     %start from middle
%     nq = size(prob_set,2)-1;
%     probability_obj_set = prob_set;
%     partworths = partworths_set;
%     pairs = pairs_set;
%     for i = 1:nq
%         query_id = queryID(min(pairs(i,:)),max(pairs(i,:)));
%         queryID(min(pairs(i,:)),max(pairs(i,:)))=0;
%         queryList(queryID==query_id)=0;
%     end
%     nquery = sum(sum(queryID>0));
    
    % find the next query
    while  nq>=1
        C = nq; % update C
        
        % calculate A
        [A,W,I,unique_I,w0] = appObjDistribution(s,d,W0,Xf,dX,dXID,c,alg,...
            pairs(1:nq-1,:),C);
        probability_obj = A/sum(A);
        fprintf('%d, %.2f, %.2f %.2f \n',nq,JSdivergence(probability_obj,target),...
            probability_obj(target_best), corr(w0',w));
        
        if nq>MAX_ITER
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
            
            actual_inq = length(aid);
            if actual_inq>0
                query_id = aid;
                options = a;
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
                CC = rho.*log2(rho)-sum(probability_obj*ones(1,actual_inq).*rho_k.*log2(rho_k))';
                
                if sum(CC.^2)<1e-16 % all C = 0
                    A = dX(DX(sub2ind(size(DX),pairs(1:nq-1,1),pairs(1:nq-1,2))),:);
                    A = bsxfun(@times,A,sign(pairs(1:nq-1,2)-pairs(1:nq-1,1)));
                    model = train(ones(size(A,1),1),...
                        sparse(A), sprintf('-s 0 -e 1e-6 -q -c %f', C));
                    w0 = (model.w)';
%                     w0 = (A'*A+1/C*eye(size(A,2)))\A'*ones(size(A,1),1);
                    B = (eye(size(A,2))-w0*w0'/(w0'*w0))*(A'*A+1/C*eye(size(A,2)));
                    [V,D] = eig(B);
                    e = diag(D);
                    v = V(:,e==min(e(e>1e-3)));
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
                    w0 = w0';
                else
                    query_id = aid(find(CC==min(CC),1));
                    options = a(find(CC==min(CC),1),:);
                end

                if isempty(options)
                    probability_obj_set(:,nq) = probability_obj;
                    partworths(:,nq) = w0';
                    fprintf('no feasible query. terminated.');
                    return;
                end
                queryID(min(options),max(options))=0;
                queryList(queryID==query_id)=0;
                
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
                nquery = nquery - 1;
            else
                probability_obj_set(:,nq) = probability_obj;
                partworths(:,nq) = w0;
                fprintf('no feasible query. terminated.');
                return;
            end
        end
    end