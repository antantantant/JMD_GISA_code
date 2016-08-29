function C = queryObj(W,I,unique_I,dx,nx,probability_obj)
    temp_set = zeros(nx,1);
    u = W*dx>0; % get all utility signs
    u = bsxfun(@times, u, I); % label all positive signs with query numbers
    for j = 1:length(unique_I)
        temp_set(unique_I(j),:) = sum(u==unique_I(j),1)/sum(I==unique_I(j));
    end

    probability_query_pos = sum(bsxfun(@times,temp_set,probability_obj),1)';% probability of new query being positive

    rho_k = zeros(nx,1);
    rho = zeros(1,1);
    probability_query_pos(probability_query_pos==0)=1e-99;

    rho_k(temp_set>0.5) = temp_set(temp_set>0.5);
    rho_k(temp_set<=0.5) = 1-temp_set(temp_set<=0.5);
    rho(probability_query_pos>0.5) = probability_query_pos(probability_query_pos>0.5);
    rho(probability_query_pos<=0.5) = 1-probability_query_pos(probability_query_pos<=0.5);
    C = rho.*log2(rho)-sum(probability_obj.*rho_k.*log2(rho_k))';