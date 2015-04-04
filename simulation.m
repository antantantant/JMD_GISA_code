clear;
close all;
addpath('..\..\..\Code\Tools\liblinear\matlab');
nexp = 100;
% design alternatives
%size 1 2 3, mem 1 2 3, keyboard 1 2, cpu 1 2, battery 1 2, price 1
%size 7 10 13, mem 32 64 128, keyboard no yes, cpu low high, battery half
%full
X = [0 1 0 1 0 0 1 0 1 0 0 1 600;
     0 1 0 1 0 0 0 1 1 0 0 1 700;
     0 1 0 0 1 0 1 0 1 0 0 1 700;
     0 1 0 0 1 0 0 1 1 0 0 1 800;
     0 1 0 0 0 1 1 0 1 0 0 1 800;
     0 1 0 0 0 1 0 1 1 0 0 1 900;
     0 0 1 1 0 0 1 0 1 0 0 1 750;
     0 0 1 1 0 0 0 1 1 0 0 1 850;
     0 0 1 0 1 0 1 0 1 0 0 1 850;
     0 0 1 0 1 0 0 1 1 0 0 1 950;
     0 0 1 0 0 1 1 0 1 0 0 1 950;
     0 0 1 0 0 1 0 1 1 0 0 1 1050;
     0 1 0 1 0 0 1 0 0 1 0 1 900;
     0 1 0 1 0 0 0 1 0 1 0 1 1000;
     0 1 0 0 1 0 1 0 0 1 0 1 1000;
     0 1 0 0 1 0 0 1 0 1 0 1 1100;
     0 1 0 0 0 1 1 0 0 1 0 1 1100;
     0 1 0 0 0 1 0 1 0 1 0 1 1200;
     0 0 1 1 0 0 1 0 0 1 0 1 1050;
     0 0 1 1 0 0 0 1 0 1 0 1 1150;
     0 0 1 0 1 0 1 0 0 1 0 1 1150;
     0 0 1 0 1 0 0 1 0 1 0 1 1250;
     0 0 1 0 0 1 1 0 0 1 0 1 1250;
     0 0 1 0 0 1 0 1 0 1 0 1 1350;
     0 1 0 1 0 0 1 0 1 0 1 0 450;
     0 1 0 1 0 0 0 1 1 0 1 0 550;
     0 1 0 0 1 0 1 0 1 0 1 0 550;
     0 1 0 0 1 0 0 1 1 0 1 0 650;
     0 1 0 0 0 1 1 0 1 0 1 0 650;
     0 1 0 0 0 1 0 1 1 0 1 0 750;
     0 0 1 1 0 0 1 0 1 0 1 0 600;
     0 0 1 1 0 0 0 1 1 0 1 0 700;
     0 0 1 0 1 0 1 0 1 0 1 0 700;
     0 0 1 0 1 0 0 1 1 0 1 0 800;
     0 0 1 0 0 1 1 0 1 0 1 0 800;
     0 0 1 0 0 1 0 1 1 0 1 0 900;
     0 1 0 1 0 0 1 0 0 1 1 0 750;
     0 1 0 1 0 0 0 1 0 1 1 0 850;
     0 1 0 0 1 0 1 0 0 1 1 0 850;
     0 1 0 0 1 0 0 1 0 1 1 0 950;
     0 1 0 0 0 1 1 0 0 1 1 0 950;
     0 1 0 0 0 1 0 1 0 1 1 0 1050;
     0 0 1 1 0 0 1 0 0 1 1 0 900;
     0 0 1 1 0 0 0 1 0 1 1 0 1000;
     0 0 1 0 1 0 1 0 0 1 1 0 1000;
     0 0 1 0 1 0 0 1 0 1 1 0 1100;
     0 0 1 0 0 1 1 0 0 1 1 0 1100;
     0 0 1 0 0 1 0 1 0 1 1 0 1200;];
% X = [X,bsxfun(@times,X(:,1:12),X(:,end))];
X(:,[1,2,7,9,11])=[];%remove 7inch, all binary ones only need one variable
Xstd = max(X(:,end));
X(:,end) = bsxfun(@rdivide, X(:,end), max(X(:,end)));
Xmean = mean(X(:,end));
X(:,end) = bsxfun(@minus, X(:,end), mean(X(:,end)));
X = [X,bsxfun(@times,X(:,1),X(:,5))];
X = [X,bsxfun(@times,X(:,1),X(:,6))];
X = [X,bsxfun(@times,X(:,1),X(:,7))];
% X = projectSphere(X);
% X(:,13) = (X(:,13)-mean(X(:,13)))/std(X(:,13))/2;
nx = size(X,1);
d = size(X,2);
cost_appGGBS = zeros(1,nexp);
cost_EGO = zeros(1,nexp);

for iter = 1:nexp
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
%% approximated GGBS
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
    model = train(ones(size(Z,1),1),Z,'-s 0 -c 1e6 -e 1e-6 -q');
    wh = model.w';
%     wh = Fan(Z);
    wh = wh/norm(wh);
    infoappGGBS.theta(i)=max(min(Z*wh),0)^d;
end
infoappGGBS.theta = infoappGGBS.theta/sum(infoappGGBS.theta);
infoappGGBS.area = infoappGGBS.theta;
infoappGGBS.nq = nq;
infoappGGBS.candidate = 1:infoappGGBS.nx;
[E_appGGBS,wh,pairs_appGGBS] = appGGBS(...
    infoappGGBS,[],w(1,:)',bestID,0);

%% EGO
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
cost_appGGBS(iter) = size(pairs_appGGBS,1);
cost_EGO(iter) = size(pairs_EGO,1);
end