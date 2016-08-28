function [f,W,I,unique_I,w0,C,weights,expected_value] = appObjDistribution(s,d,W0,X,c,A,competitors)
    [W,w0,C,weights] = sampling(s,d,A);
    
    util = (W*X');
    num_competitor = length(competitors);
    util_competitor = util(:,competitors);
    util_competitor_all = kron(util_competitor,ones(1,size(X,1)));
    util_all = repmat(util,1,num_competitor);
    exp_delta_util = exp(-util_all+util_competitor_all);
    exp_sum_delta_util = exp_delta_util*sparse(repmat(eye(size(X,1)),num_competitor,1))+1;
    obj_app = bsxfun(@plus,-log(exp_sum_delta_util),log(c'));
    
%     obj_app = bsxfun(@plus,-log(1+exp(-util)),log(c'));
    
    best = bsxfun(@eq, obj_app, max(obj_app,[],2));
    best = best.*util;
    best(best==0) = -1e32;
    best = bsxfun(@eq, best, max(best,[],2));
    f = (weights'*best/sum(weights))';

    % each line of best needs to have only one 1
    if any(sum(best,2)>1)
        error = 1;
    end

    I = best*(1:size(best,2))';
    %     obj = bsxfun(@times,1./(1+exp(-W*X(keep,:)')),c(keep)');
    %     f(keep) = sum(bsxfun(@eq, obj, max(obj,[],2))/s)';
        
    unique_I = unique(I);
    unique_I(unique_I>length(c))=[];
    
    % to calculate expected profit of each candidate
    expected_value = weights'*exp(obj_app)/sum(weights);
    