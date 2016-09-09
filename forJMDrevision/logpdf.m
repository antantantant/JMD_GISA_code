function p = logpdf(x,A,C)
    dW = sum(x.^2,2);
    temp = -x*A';
    temp2 = log(1+exp(temp));
    temp2(temp2==inf) = temp(temp2==inf);
    p = -dW/2/C+sum(-temp2,2);