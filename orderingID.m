function id = orderingID(info,query,query_id)
    check = info.query_index(:,query_id)*diag(query);
    id = sum(check,2)==length(query);
%     
%     for j = 1:info.nbeta
%         o = info.beta_unique_label(j,:);
%         flag = true;
%         for i = 1:size(pair,1)
%             if sum(o(find(o==pair(1)):end)==pair(2))==0
%                 flag = false;
%                 break;
%             end
%         end
%         if flag
%             id(j) = true;
%         end
%     end