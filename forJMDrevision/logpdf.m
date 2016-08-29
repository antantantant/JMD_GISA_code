function f = logpdf(x,A,C,theta)
    f = -sum(log(1+exp(-theta*A*x))) - 1/2/C*(x'*x);