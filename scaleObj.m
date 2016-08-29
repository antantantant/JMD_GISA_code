function f = scaleObj(x,y,w,z)
    cv = 3;
    ci = 1e6;
    s = 5e6;
    theta = 1;
    a = 1/(1+exp(-theta*calMarketShare(x,y,w,z)));
    if isempty(x)
        p = z(6);
    else
        p = x(15);
    end
    f = -(s*a*(p-cv)-ci);