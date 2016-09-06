% output optimal C and associated w
function [w, best_C] = crossvalidation(X,y)

%%%%%%%%%%%%%%%% ALL C, NO CROSSVALIDATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     fprintf('\n svm training...');
%     [n,p] = size(X);
%     
%     A = bsxfun(@times,X,y);
%     C = 1/n*[1,5,1e1,5e1,1e2,5e2,1e3,5e3,1e4,5e4,1e5,5e5,1e6,5e6,1e7,5e7,1e8];
%     weights = zeros(1,length(C));
%     w = zeros(length(C),p);
%     for i = 1:length(C)
%         model = train(sparse(ones(n,1)), sparse(A), ...
%             sprintf('-s 0 -e 1e-3 -q -c %f', C(i)));
%         w(i,:) = model.w;
%         weights(i) = likelihood(model.w,C(i),A);
%     end
%     weights = weights/sum(weights);
%     
%     function f = likelihood(w,C,A)
%         f = exp(-1/2/C*(w*w'))*prod(1./(exp(-A*w')+1));
    
%     best_C = C;
%     weights = cv_acc/sum(cv_acc);
%     w = zeros(length(best_C),p);
%     for i = 1:length(best_C)
%         model = train(sparse([ones(n,1);-ones(n,1)]), sparse([A;-A]), ...
%             sprintf('-s 0 -e 1e-3 -q -c %f', best_C(i)));
%         w(i,:) = model.w;
%     end
    
    
%     %# pair (C,gamma) with best accuracy
%     idx = find(cv_acc==max(cv_acc));
%     %# now you can train you model using best_C and best_gamma
%     best_Cs = C(idx);
%     best_C = unique([min(best_Cs), max(best_Cs)]);
%     
%     w = zeros(length(best_C),p);
%     for i = 1:length(best_C)
%         model = train(sparse([ones(n,1);-ones(n,1)]), sparse([A;-A]), ...
%             sprintf('-s 0 -e 1e-3 -q -c %f', best_C(i)));
%         w(i,:) = model.w;
%     end
%     fprintf('done. \n');


%%%%%%%%%%%%%%%% CROSSVALIDATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     fprintf('\n svm training...');
    [n,p] = size(X);
    
    
    folds = 10;
    C = [1e-1,1,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8]; % search for C
%     C = [10000];
    
    %# grid search, and cross-validation
    cv_acc = zeros(numel(C),1);
    A = bsxfun(@times,X,y);
    for i=1:numel(C)
%         cv_acc(i) = train(sparse([ones(n,1);-ones(n,1)]), sparse([A;-A]), ...
        cv_acc(i) = train(sparse(ones(n,1)), sparse(A), ...
        sprintf('-s 0 -e 1e-6 -q -c %f -v %d', C(i), folds));
    end
    
    %# pair (C,gamma) with best accuracy
    idx = find(cv_acc==max(cv_acc));
    idx = idx(end);
    
    %# now you can train you model using best_C and best_gamma
    best_C = C(idx);
    
%     model = train(sparse([ones(n,1);-ones(n,1)]), sparse([A;-A]), ...
    model = train(sparse(ones(n,1)), sparse(A), ...
        sprintf('-s 0 -e 1e-6 -q -c %f', best_C));
    w = (model.w);
%     w = ((A'*A+1/best_C*eye(size(A,2)))\A'*ones(size(A,1),1))';
%     weights = 1;
%     fprintf('done. \n');

%     for i = 1:length(C)
%         w = ((A'*A+1/C(i)*eye(size(A,2)))\A'*ones(size(A,1),1))';
%         f = likelihood(w,C(i),A)
%     end
    
    function f = likelihood(w,C,A)
        f = exp(-1/2/C*(w*w'))*prod(exp(-sum((A*w'-1).^2)));
    