TEST = 12;
target_best_set = zeros(TEST,1);
for test = 1:TEST
    rng(test);
    fprintf('\n%%%%%%%% test number %d %%%%%%%%%%',test);
    
    % Calculate true distribution under sigma
    Wtrue = mvnrnd(wtrue,Dw,s);
    util = (Wtrue*Xf(1:nt,:)');
    num_competitior = 10;
    competitors = randperm(size(Xf,1),num_competitior);
    util_competitor = util(:,competitors);
    util_competitor_all = kron(util_competitor,ones(1,size(Xf,1)));
    util_all = repmat(util,1,num_competitior);
    exp_delta_util = exp(-util_all+util_competitor_all);
    exp_sum_delta_util = exp_delta_util*repmat(eye(size(Xf,1)),num_competitior,1);
    obj_app = bsxfun(@plus,-log(exp_sum_delta_util),log(c(1:nt,:)'));
    best = bsxfun(@eq, obj_app, max(obj_app,[],2));
    best = best.*util;
    best(best==0) = -1e9;
    best = bsxfun(@eq, best, max(best,[],2));
    target_dist = sum(best/s)';
    [~,target_best] = max(target_dist);
    target_best_set(test) = target_best;
end
