function [w_id,unique_w_id,w_ii] = getXstar(w,X,r)
    w_id = zeros(size(w,1),r);
    w_ii = zeros(size(w,1),1);
    for i = 1:size(w,1)
        [u,uid] = sort(w(i,:)*X','descend');
        w_id(i,:) = uid(1:r);
    end
    unique_w_id = unique(w_id,'rows');
    for i = 1:size(w,1)
        for j = 1:size(unique_w_id)
            if sum(abs(w_id(i,:)-unique_w_id(j,:)))==0
                w_ii(i) = j;
                break;
            end
        end
    end