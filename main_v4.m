clear;
close all;
matlabpool open;
addpath('..\..\..\Code\Tools\liblinear\matlab');
nexp = 100;
cost_GGBS = zeros(10,nexp);
cost_appGGBS = zeros(10,nexp);
cost_EGO = zeros(10,nexp);
time_GGBS = zeros(10,nexp);
time_appGGBS = zeros(10,nexp);
time_EGO = zeros(10,nexp);

for nx = 10
    for iter = 1:nexp
    % design alternatives
    d = 5; 
    X = randn(nx,d);
    X = projectSphere(X);
    % 
    w = rand(1,d); % try different distributions
    w = projectSphere(w);
    bestID = find(w*X'==max(w*X'));
    % load test.mat;
    % comparison matrix
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
    % learning parameters
    par.C = 1;

    %% GGBS
%     t = cputime;
%     sample_size = 1e3;
%     sample = rand(sample_size,d)*2-1;
%     sample = projectSphere(sample);
%     [temp,sample_label] = sort(sample*X',2,'descend');
%     unique_label = unique(sample_label,'rows');
%     no = size(unique_label,1);
%     sample_label_id = zeros(sample_size,1);
%     beta = zeros(no,1);
%     for i = 1:sample_size
%         disp(num2str(i));
%         for j = 1:size(unique_label)
%             if sum(abs(sample_label(i,:)-unique_label(j,:)))==0
%                 sample_label_id(i) = j;
%                 beta(j) = beta(j) + 1;
%                 break;
%             end
%         end
%     end
%     beta = beta/sum(beta);
%     index = -ones(no,size(dX,1));
%     for j = 1:no
%         disp(num2str(j));
%         o = unique_label(j,:);
%         for i = 1:size(dX,1)
%             [obj1,obj2] = find(dXID==i);
%             if sum(o(find(o==obj1):end)==obj2)>0
%                 index(j,i) = 1;
%             end
%         end
%     end
%     infoGGBS.X = X;
%     infoGGBS.nx = nx;
%     infoGGBS.d = d;
%     infoGGBS.r = 1;
%     infoGGBS.beta = beta;
%     infoGGBS.nbeta = no;
%     infoGGBS.beta_unique_label = unique_label;
%     infoGGBS.query_index = index;
%     infoGGBS.Z = dX;
%     infoGGBS.Z_id = dXID;
%     infoGGBS.sample_size = sample_size;
%     infoGGBS.sample_label = sample_label(:,1);
%     infoGGBS.unique_label = (1:nx)';
%     infoGGBS.no = size(infoGGBS.unique_label,1);
%     infoGGBS.sample_label_id = infoGGBS.beta_unique_label(:,1);
%     infoGGBS.theta = zeros(infoGGBS.no,1);
%     for j = 1:nx
%         infoGGBS.theta(j) = sum(infoGGBS.sample_label==j);
%     end
%     infoGGBS.scale = sum(infoGGBS.theta);
%     infoGGBS.area = infoGGBS.theta;
%     infoGGBS.theta = infoGGBS.theta/sum(infoGGBS.theta);
%     [E_GGBS,wh,pairs_GGBS] = GGBS(infoGGBS,[],[],w',bestID,0);
%     if isempty(wh) 
%         continue;
%     end
%     time_GGBS(nx/5,iter) = cputime - t;

    %% approximated GGBS
    t = cputime;
    infoappGGBS.X = X;
    infoappGGBS.dX = dX;
    infoappGGBS.dXID = dXID;
    infoappGGBS.nx = nx;
    infoappGGBS.d = d;
    infoappGGBS.theta = zeros(infoappGGBS.nx,1);
    for i = 1:nx
        Z = zeros(nx-1,d);
        count = 1;
        for j = [1:i-1,i+1:nx]
            Z(count,:) = X(i,:)-X(j,:);
            count = count + 1;
        end
        Z = sparse(Z);
        model = train(ones(size(Z,1),1),Z,'-s 0 -c 1e6 -e 1e-6');
        wh = model.w';
        wh = wh/norm(wh);
        infoappGGBS.theta(i)=max(min(Z*wh),0)^2;
    end
    infoappGGBS.scale = sum(infoappGGBS.theta);
    infoappGGBS.theta = infoappGGBS.theta/sum(infoappGGBS.theta);
    infoappGGBS.nq = nq;
    infoappGGBS.candidate = 1:infoappGGBS.nx;
    [E_appGGBS,wh,pairs_appGGBS] = appGGBS(...
        infoappGGBS,[],w(1,:)',bestID,0);
    time_appGGBS(nx/5,iter) = cputime - t;

    %% EGO
    t = cputime;
    infoEGO.X = X;
    infoEGO.dX = dX;
    infoEGO.dXID = dXID;
    infoEGO.nx = nx;
    infoEGO.d = d;
    infoEGO.theta = ones(infoEGO.nx,1)/infoEGO.nx;
    infoEGO.nq = nq;
    infoEGO.candidate = 1:infoEGO.nx;
    [E_EGO,wh,pairs_EGO] = EGO(...
        infoEGO,[],w(1,:)',bestID,0);
    time_EGO(nx/5,iter) = cputime - t;

    %%
    cost_GGBS(nx/5,iter) = size(pairs_GGBS,1);
    cost_appGGBS(nx/5,iter) = size(pairs_appGGBS,1);
    cost_EGO(nx/5,iter) = size(pairs_EGO,1);
    end
%     saveto = ['nx',num2str(nx),'_d',num2str(d),'.mat'];
%     save(saveto);
end
mean(cost_GGBS)
mean(cost_appGGBS)
mean(cost_EGO)
