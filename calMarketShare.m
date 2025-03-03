function f = calMarketShare(x,y,w,z)
    if isempty(z)
        z1 = 4*pi*x(6)*x(9)*x(10)*(x(1)+x(2))*(x(3)+x(4))...
            /x(11)/(x(1)*(x(3)+x(4))+x(3)*(x(1)+x(5)));
        z2 = x(13)/x(14);
        z3 = x(13)*x(14);
        z4 = pi*x(12)/z1;
        z5 = (2*tan(pi*y(11)/z1))*(x(12)/2-y(10))/(1+2/y(12)*tan(pi*(y(11)/z1)));
        p = x(15);
    else
        z1 = z(1);
        z2 = z(2);
        z3 = z(3);
        z4 = z(4);
        z5 = z(5);
        p = z(6);
    end
    
    W1 = w(1:5);
    W2 = w(6:10);
    W3 = w(11:15);
    W4 = w(16:20);
    W5 = w(21:25);
    W6 = w(26:30);
    
    Z1 = 200:50:400;
    Z2 = [0.75,0.88,1.00,1.14,1.33];
    Z3 = 100:10:140;
    Z4 = 1/16:1/32:6/32;
    Z5 = 0.75:0.25:1.75;
    Z6 = 10:5:30;
    
    w1 = spline(Z1,W1,z1);
    w2 = spline(Z2,W2,z2);
    w3 = spline(Z3,W3,z3);
    w4 = spline(Z4,W4,z4);
    w5 = spline(Z5,W5,z5);
    w6 = spline(Z6,W6,p);
    
    f = w1+w2+w3+w4+w5+w6;
   