function [Pr_set,Pq] = updatePr(info,q_id,q,nx)
    Pr_set = zeros(nx,1);
    if ~isempty(q_id)
        Pq = zeros(info.No,1); % calculate the probability of the current query response
        for j = 1:info.No
            q_o = info.q_o{j};
            if isempty(q_o)
                Pq(j) = 0;
            else
                q_o = q_o(:,q_id);
                t = q_o - ones(size(q_o,1),1)*q;
                Pq(j) = info.Po(j)*sum(diag(t*t')==0)/size(q_o,1);
            end
        end
        Pr_set = Pq/sum(Pq);
    else
        Pr_set = info.Po;
        Pq = [];
    end