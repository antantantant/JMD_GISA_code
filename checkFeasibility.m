function [design_set,feasible] = checkFeasibility(X,Z)
    
    y1 = 0.3;
    y2 = 0.5;
    y3 = 1.9;
    y4 = 0.25;
    y5 = 3;
    y6 = 2;
    y7 = 1.13;
    y8 = 1;
    y9 = 1.1;
    y10 = 0.31;
    y11 = 16;
    y12 = 1.29;
    y13 = 4;
    y = [y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13];
    % linear constraints
    A = zeros(6,14);
    b = zeros(6,1);
    A(1,[12,14]) = [1,-1];
    b(1) = -2*y1;
    A(2,[7,12,13]) = [1,1,-1];
    b(2) = -2*y1 - y9;
    A(3,[4,5,13]) = [1,1,-1];
    b(3) = -2*y1;
    A(4,[2,5]) = [-1,1];
    A(5,[7,8,11,13]) = [1,1,1,-1];
    b(5) = -2*y1 - y9;
    A(6,[7,8,10,12,13]) = [1,1,1,0.5,-1];
    b(6) = -2*y1 - y7 - y9;
    % solve
    lb = [0.125*ones(5,1);1;0.5;1;0.25;0.5;0.5;1;1;1];
    ub = [36;36;24;24;36;200;12;36;24;1.9;1.9;36;36;36];
    load trueSolution.mat;
    x = x(1:14);
    X = X(:,1:25);
    [X,IA,IC] = unique(X,'rows');
    feasible = ones(size(X,1),1);
    design_set = nan(size(X,1),length(x));
    for i = 1:size(X,1)
        disp(num2str(i));
        bb = X(i,:);
        z = Z(bb==1);
        [xx,fval,exitflag] = fmincon(@(x)eqn(x,z,y),x,A,b,[],[],lb,ub,@(x)scaleCon(x,y));
%         options = gaoptimset('InitialPopulation',x);
%         [x,fval,exitflag] = ga(@(x)eqn(x,z,y),14,A,b,[],[],lb,ub,@(x)scaleCon(x,y),options)
        if (fval>1e-3 && exitflag>0) || exitflag<0
            feasible(i)=0;
        end
        design_set(i,:) = xx';
    end
    feasible = feasible(IC);
    design_set = design_set(IC,:);
    disp('finished feasibility check.');
    
    