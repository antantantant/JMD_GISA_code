function x_new = proprnd(x)
    x_new = x + (rand(1,size(x,2)) - 0.5);  