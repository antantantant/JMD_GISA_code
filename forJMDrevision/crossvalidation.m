% output optimal C and associated w
function [w, best_C] = crossvalidation(X,Y)
    fprintf('\n svm training...');
    [n,p] = size(X);
    
    folds = 5; % 5-fold cv
    C = 1/n*[1e1,1e2,1e3]; % search for C
    %# grid search, and cross-validation
    cv_acc = zeros(numel(C),1);
    for i=1:numel(C)
        cv_acc(i) = train(Y, X, ...
        sprintf('-s 0 -e 1e-9 -q -c %f -v %d', C(i), folds));
    end
    %# pair (C,gamma) with best accuracy
    [~,idx] = max(cv_acc);
    %# now you can train you model using best_C and best_gamma
    best_C = C(idx);
    
    model = train(Y, X, ...
        sprintf('-s 0 -e 1e-6 -q -c %f', best_C));
    w = model.w;
    fprintf('done. \n');
    