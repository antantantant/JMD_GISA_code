% clear;
% close all;
addpath('..\..\..\Code\Tools\liblinear\matlab');

nexp = 1;
d = 10;
nx = 100;
nn = 2;

for exp = 1:nexp
X = randn(nx,d);
X = projectSphere(X);
w = randn(nn,d); % try different distributions
w = projectSphere(w);
bestID = zeros(nn,1);
for i = 1:nn
    bestID(i) = find(w(i,:)*X'==max(w(i,:)*X'));
end

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
    Z = bsxfun(@rdivide,Z,sqrt(dot(Z,Z,2)));
    Z = sparse(Z);
    model = train(ones(size(Z,1),1),Z,'-s 0 -c 1e6 -e 1e-6 -q');
    wh = model.w';
%     wh = Fan(Z);
    wh = wh/norm(wh);
    infoappGGBS.theta(i)=max(min(Z*wh),0)^(d-1);
%     infoappGGBS.theta(i)=max(min(Z*wh),0);
end
infoappGGBS.scale = sum(infoappGGBS.theta);
infoappGGBS.theta = infoappGGBS.theta/sum(infoappGGBS.theta);
infoappGGBS.area = infoappGGBS.theta;
infoappGGBS.nq = nq;
infoappGGBS.candidate = 1:infoappGGBS.nx;
infoEGO.X = X;
infoEGO.dX = dX;
infoEGO.dXID = dXID;
infoEGO.nx = nx;
infoEGO.d = d;
infoEGO.scale = infoappGGBS.scale;
infoEGO.theta = infoappGGBS.theta;
infoEGO.area = infoappGGBS.area;
infoEGO.nq = nq;
infoEGO.candidate = 1:infoEGO.nx;

for iter = 1:nn*50
%% approximated GGBS
[E_appGGBS,wh,pairs_appGGBS] = appGGBS(...
    infoappGGBS,[],w(mod(iter-1,nn)+1,:)',bestID(mod(iter-1,nn)+1),0);
%% EGO
[E_EGO,wh,pairs_EGO] = EGO(...
    infoEGO,[],w(mod(iter-1,nn)+1,:)',bestID(mod(iter-1,nn)+1),0);
cost_appGGBS(iter) = size(pairs_appGGBS,1);
cost_EGO(iter) = size(pairs_EGO,1);
infoappGGBS = updateInfo(infoappGGBS,bestID(mod(iter-1,nn)+1),iter);
infoEGO = updateInfo(infoEGO,bestID(mod(iter-1,nn)+1),iter);
end
saveto = ['exp',num2str(exp),'_update.mat'];
save(saveto);
end