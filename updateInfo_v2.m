function info = updateInfo_v2(info)
    w = info.wh; %learned w from observations
    
    % for histogram
    [t,id] = sort(w(end,:)*info.X',2,'descend');
    theta_n = info.theta*info.norm;
    theta_n(id(1)) = theta_n(id(1)) + info.tau;
    info.theta = theta_n/sum(theta_n);
    info.theta_set(:,size(w,1)+1) = info.theta;