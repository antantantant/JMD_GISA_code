function [f,W] = appObjDistribution(s,d,W,Xf,dX,dXID,c,const,alg,pairs,p)
    nx = size(Xf,1);
    X = Xf;
    f = zeros(nx,1);
    
    if (strcmp(alg,'pref_only'))
%         % appggbs for preference only
%         % calculate areas
%         f = appDistribution(nx,Xf,d,pairs,p);
%         % calculate volume
%         f = f.*const;
        if ~isempty(p)
            keep = p>1e-3/nx;
        else
            keep = true(nx,1);
        end
        dXID = dXID + dXID';

        W = sampling(s,d,W,Xf,dX,dXID,pairs);
        try
            obj = W*X(keep,:)';
        catch
            wait = 1;
        end
        f(keep) = sum(bsxfun(@eq, obj, max(obj,[],2))/s)';
        
    else
        % proposed sampling method
        if ~isempty(p)
            keep = p>1e-3/nx;
        else
            keep = true(nx,1);
        end
        dXID = dXID + dXID';

        W = sampling(s,d,W,Xf,dX,dXID,pairs);
        
%         if sum(feasible)>1e3
%             segl = 1e3;
%             segn = ceil(sum(feasible)/segl);
%             for i = 1:segn
%                 id = ((i-1)*segl+1):min(i*segl,sum(feasible));
%                 obj = bsxfun(@times,1./(1+exp(-Wt(id,:)*X(keep,:)')),c(keep)');
%                 f(keep) = f(keep) +...
%                     sum(bsxfun(@times, bsxfun(@eq, obj, max(obj,[],2))/s, ht(id)))';
%             end
% 
%         else
            obj = bsxfun(@times,1./(1+exp(-W*X(keep,:)')),c(keep)');
%             f(keep) = sum(bsxfun(@times, bsxfun(@eq, obj, max(obj,[],2))/s, h))';
            f(keep) = sum(bsxfun(@eq, obj, max(obj,[],2))/s)';

%         end
    end
    