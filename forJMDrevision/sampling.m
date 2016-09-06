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
        W = mvnrnd(zeros(1,d),Sigma,s);
        qW = log(mvnpdf(W, zeros(1,d), Sigma));
        dW = sum(W.^2,2);
        pW = -dW/2/C+sum(-log(1+exp(-W*A')),2);
        weights = exp(pW-qW);
%         temp = -W*A';
%         temp2 = log(1+exp(temp));
%         temp2(temp2==inf) = temp(temp2==inf);
%         pW = -dW/2/C+sum(-temp2,2);
%         weights = exp(pW-qW);
        if sum(weights)==0
            weights=ones(s,1);
        end
        weights = weights/sum(weights);
    end