function info = updateInfo(info,bestID,iter)
    
    info.theta = info.theta*(info.nx+(iter-1)*1);
    info.theta(bestID) = 1+info.theta(bestID);
    info.theta = info.theta/(info.nx+iter*1);
    
