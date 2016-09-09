function y = bisectfunc(x, w0, Sigma, A, C)
%     qW = log(mvnpdf(w0, w0, Sigma*x));
%     dW = sum(w0.^2,2);
%     temp = -w0*A';
%     temp2 = log(1+exp(temp));
%     temp2(temp2==inf) = temp(temp2==inf);
%     pW = -dW/2/C+sum(-temp2,2);
%     pqW = pW-qW-log(x);
    s = 1e4;
    W = mvnrnd(w0,Sigma*x,s);
    qW = log(mvnpdf(W, w0, Sigma*x));
    dW = sum(W.^2,2);
    temp = -W*A';
    temp2 = log(1+exp(temp));
    temp2(temp2==inf) = temp(temp2==inf);
    pW = -dW/2/C+sum(-temp2,2)-log(x);
    pqW = pW-qW;
    a = max(pqW);
    y = a+log(sum(exp(bsxfun(@minus,pqW,a))))-log(s);