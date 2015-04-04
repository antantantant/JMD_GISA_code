function [f,W,I,unique_I,w0] = appObjDistribution(s,d,W0,Xf,dX,dXID,c,alg,pairs,C)
    nx = size(Xf,1);
    X = Xf;
    f = zeros(nx,1);

    dXID = dXID + dXID';

    [W,w0] = sampling(s,d,W0,Xf,dX,dXID,pairs,C);
    
    util = (W*X');
%     if strcmp(alg,'pref_only')
%         obj_app = util;
%     else
        obj_app = bsxfun(@plus,-log(1+exp(-util)),log(c'));
%     end
        
    best = bsxfun(@eq, obj_app, max(obj_app,[],2));
    best = best.*util;
    best(best==0) = -1e9;
    best = bsxfun(@eq, best, max(best,[],2));
    f = sum(best/s)';
    
    % each line of best needs to have only one 1
    if any(sum(best,2)>1)
        error = 1;
    end
    
    I = best*(1:size(best,2))';
    unique_I = unique(I);
%     obj = bsxfun(@times,1./(1+exp(-W*X(keep,:)')),c(keep)');
%     f(keep) = sum(bsxfun(@eq, obj, max(obj,[],2))/s)';
    