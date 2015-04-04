function [probability_obj_set,pairs] = ...
    bestProb(s,d,Xf,dX,dXID,c,inq,const,alg,queryID,pairs,w,target,nq)
    
    nx = size(Xf,1);
    
    % calculate A
    [A,W] = appObjDistribution(s,d,w,Xf,dX,dXID,c,const,alg,pairs,[]);
    probability_obj = A/sum(A);
    fprintf('|%d, %.2f|',nq,JSdivergence(probability_obj,target));    

    % find the next query
    if  nq>=0
        if nq>100
            probability_obj_set = probability_obj;
            return;
        else
            [~,sort_obj] = sort(probability_obj,'descend');
            if ~isempty(pairs)
                count = 1;
                sampled = unique(pairs);
                new_sample = [];
%                 while count<=inq
%                     if (sum(sampled==sort_obj(count))>0 ||...
%                             probability_obj(sort_obj(count))<(1/nx*1e-3))
%                         count = count + 1;
%                     else
%                         new_sample = sort_obj(count);
%                         break;
%                     end
%                 end
                sort_obj = setdiff(sort_obj,sampled);
                
                if isempty(sort_obj)
                    probability_obj_set = probability_obj;
                    return;
                else
                    new_sample = sort_obj(1);
                end

                options = [pairs(end,1),new_sample];
            else
                options = [sort_obj(1),sort_obj(2)];
            end
            queryID(min(options),max(options))=0;
            u = (Xf(options(1),:)-Xf(options(2),:))*w;
            if u>0
                pairs = [pairs;options(1),options(2)];
            elseif u<0
                pairs = [pairs;options(2),options(1)];
            end
            [probability_obj_sub_set,pairs] = ...
                bestProb(s,d,Xf,dX,dXID,c,inq,const,alg,queryID,pairs,w,target,nq+1);
            probability_obj_set = [probability_obj, probability_obj_sub_set];
        end
    end 