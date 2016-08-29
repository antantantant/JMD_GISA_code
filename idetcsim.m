clear;
close all;
addpath('..\..\..\Code\Tools\liblinear\matlab');
nexp = 100;

for dim = 5
for level = 3
    d = dim*level;
    nx = level^dim;
    X = zeros(nx,d);
    for j = 1:dim
            X(1,(j-1)*level+1) = 1;
    end
    for i = 2:nx
        count = 1;
        a = zeros(dim,1);
        n = i-1;
        while ~(count>dim||n==0)
            a(count) = mod(n,level);
            n = (n-a(count))/level;
            count = count + 1;
        end
        for j = 1:dim
            X(i,(j-1)*level+a(j)+1) = 1;
        end
    end
    X = bsxfun(@minus, X, mean(X));
    X = bsxfun(@rdivide, X, std(X));

    nq = nx*(nx-1)/2;
    dX = zeros(nq,d);
    dXID = zeros(nx,nx);
    count = 1;
    for i = 1:nx
        for j = i+1:nx
            dXID(i,j) = count;
            dX(count,:) = X(i,:) - X(j,:);
            count = count + 1;
        end
    end
    cost_appGGBS = zeros(1,nexp);
    cost_EGO = zeros(1,nexp);
    time_appGGBS = zeros(1,nexp);
    time_EGO = zeros(1,nexp);

    for iter = 1:nexp
        w = rand(1,d); % try different distributions
        w = projectSphere(w);
        bestID = find(w*X'==max(w*X'));

        %% approximated GGBS
        t = cputime;
        infoappGGBS.X = X;
        infoappGGBS.dX = dX;
        infoappGGBS.dXID = dXID;
        infoappGGBS.nx = nx;
        infoappGGBS.d = d;
%         infoappGGBS.theta = zeros(infoappGGBS.nx,1);
%         for i = 1:nx
%             Z = zeros(nx-1,d);
%             count = 1;
%             for j = [1:i-1,i+1:nx]
%                 Z(count,:) = X(i,:)-X(j,:);
%                 count = count + 1;
%             end
%             Z = bsxfun(@rdivide,Z,sqrt(dot(Z,Z,2)));
%             Z = sparse(Z);
%             model = train(ones(size(Z,1),1),Z,'-s 0 -c 1e6 -e 1e-6 -q');
%             wh = model.w';
%             wh = wh/norm(wh);
%             infoappGGBS.theta(i)=max(min(Z*wh),0)^d;
%         end
%         infoappGGBS.theta = infoappGGBS.theta/sum(infoappGGBS.theta);
        infoappGGBS.theta = ones(infoappGGBS.nx,1)/infoappGGBS.nx;
        infoappGGBS.area = infoappGGBS.theta;
        infoappGGBS.nq = nq;
        infoappGGBS.candidate = 1:infoappGGBS.nx;
        [E_appGGBS,wh,pairs_appGGBS] = appGGBS(...
            infoappGGBS,[],w(1,:)',bestID,0);
        time_appGGBS(iter) = cputime - t;

        %% EGO
        t = cputime;
        infoEGO.X = X;
        infoEGO.dX = dX;
        infoEGO.dXID = dXID;
        infoEGO.nx = nx;
        infoEGO.d = d;
        infoEGO.theta = infoappGGBS.theta;
        infoEGO.area = infoappGGBS.area;
        infoEGO.nq = nq;
        infoEGO.candidate = 1:infoEGO.nx;
        [E_EGO,wh,pairs_EGO] = EGO(...
            infoEGO,[],w(1,:)',bestID,0);
        time_EGO(iter) = cputime - t;

        cost_appGGBS(iter) = size(pairs_appGGBS,1);
        cost_EGO(iter) = size(pairs_EGO,1);
    end
    [h,p] = ttest2(cost_appGGBS,cost_EGO,0.05,'left','unequal')
    saveto = ['idetc_dim',num2str(dim),'_level',num2str(level),'.mat'];
    save(saveto);
end
end

