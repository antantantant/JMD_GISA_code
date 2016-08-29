function [W,w0,C,weights] = sampling(s,d,A)
    
    if isempty(A)
        W = mvnrnd(zeros(1,d),eye(d),s); % generate candidates
        w0 = zeros(1,d);
        C = 1;
        weights = ones(s,1);
    else
        % update estimation of w, row vector
        [w0, C] = crossvalidation(sparse(A),ones(size(A,1),1)); 
        % calculate sigma
        As = bsxfun(@times, A, sqrt(exp(A*w0'))./(1+exp(A*w0')));
        As(isnan(As))=0;
        Hessian = (eye(length(w0))/C+(As'*As));
        Sigma = inv(Hessian);
        % importance sampling
        W = mvnrnd(w0,Sigma,s);
        qW = mvnpdf(W, w0, Sigma)+1e-36;
        dW = sum(W.^2,2);
        pW = exp(-dW/2/C).*prod(1./(1+exp(-W*A')),2)+1e-36;
        weights = pW./qW;
    end