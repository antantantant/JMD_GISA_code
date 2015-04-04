% Jensen-Shannon divergence
function d = JSdivergence(a,b)
    d = sum(-0.5*(a+b).*log(0.5*(a+b+eps))+...
        -0.5*a.*log(a+eps)-0.5*b.*log(b+eps));