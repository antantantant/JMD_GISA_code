% function test_logistic
    p = 10;
    n = 1000;
    theta = 10;
    w = (rand(p,1)*2-1)*theta;
    X = rand(n,p)*2-1;
    r = rand(n,1);
    y = bsxfun(@le, r, 1./(1+exp(-X*w)))*2-1;
    error = sum(abs((y+1)/2-(X*w>0)))
    
%     opt = optimset('Algorithm', 'active-set', 'GradObj','on', 'Display', 'off');
%     w0 = fminunc(@(w)likelihood(w,X,y),zeros(p,1),opt);
    
    folds = 10;
    C = 1/n*[1e-2,1e-1,1,1e1,1e2,1e3,5e3,1e4,5e4,1e5,5e5,1e6]; % search for C
    cv_acc = zeros(numel(C),1);
    A = bsxfun(@times,X,y);
    for i=1:numel(C)
        cv_acc(i) = train(sparse([ones(n,1);-ones(n,1)]), sparse([A;-A]), ...
        sprintf('-s 0 -e 1e-6 -q -c %f -v %d', C(i), folds));
    end
    
    [~,idx] = max(cv_acc);
    best_C = C(idx);
    
    model = train(sparse([ones(n,1);-ones(n,1)]), sparse([A;-A]), ...
        sprintf('-s 0 -e 1e-6 -q -c %f', best_C));
    w0 = (model.w)';
    
    
    corr(w,w0)
    norm(w-w0)
    
%     function [f,g] = likelihood(w,X,y)
%             f = -sum(log(1+exp(-X*w.*y)));
%             g = (bsxfun(@times,X,y))'*(1./(1+exp(X*w.*y)));