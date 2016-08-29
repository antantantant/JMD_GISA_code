% % scale case study based on J. Michalek's papers
% % load liblinear
% addpath('..\..\..\Code\Tools\liblinear\matlab');
% 
% % fixed partworth from Linking Marketing and Engineering Product Design Decisions
% % via Analytical Target Cascading 
% w = [-0.534,0.129,0.228,0.104,0.052,...
%     -0.058,0.253,0.278,-0.025,-0.467,...
%     0.015,-0.098,0.049,0.047,-0.033,...
%     -0.366,-0.164,0.215,0.194,0.100,...
%     -0.744,-0.198,0.235,0.291,0.396,...
%     0.719,0.482,0.054,-0.368,-0.908];
% % attributes and levels
% z1 = 200:50:400;
% z2 = [0.75,0.88,1.00,1.14,1.33];
% z3 = 100:10:140;
% z4 = 1/16:1/32:6/32;
% z5 = 0.75:0.25:1.75;
% z6 = 10:5:30;
% z = [z1,z2,z3,z4,z5,z6];
% %%
% nq = 10; %number of queries per user
% nc = 5; %number of users
% wn = w/norm(w); %normalized w
% saveto = [num2str(nc),'nc',num2str(nq),'nq_1000candidate_deterministic.mat'];
% % if exist(saveto,'file')
% %     load(saveto);
% % end
% 
% % get w estimation from full information
% X = zeros(5^6,30);
% count = 1;
% for i1 = 1:5
%     for i2 = 1:5
%         for i3 = 1:5
%             for i4 = 1:5
%                 for i5 = 1:5
%                     for i6 = 1:5
%                         X(count,[i1,5+i2,10+i3,15+i4,20+i5,25+i6]) = 1;
%                         count = count + 1;
%                     end
%                 end
%             end
%         end
%     end
% end
% % Z = [];
% % count = 1;
% % for i = 1:(5^6-1)
% %     for j = (i+1):5^6
% %         if w*(X(i,:)-X(j,:))'>0
% %             Z(count,:) = X(i,:)-X(j,:);
% %             count = count + 1;
% %         else
% %             Z(count,:) = X(j,:)-X(i,:);
% %             count = count + 1;
% %         end
% %     end
% % end
% feasible = checkFeasibility(X,z);
% Xf = zeros(sum(feasible)*5,30);
% id = find(feasible>0);
% for i = 1:length(id)
%         Xf((i-1)*5+(1:5),:) = X((id(i)-1)*5+(1:5),:);
% end
% %%
% maxiter = 50;
% drank1 = zeros(maxiter,1); drank2 = drank1; drank3 = drank1;
% w1 = zeros(maxiter,30); w2 = w1; w3 = w1;
% rho1 = drank1; rho2 = drank1; rho3 = drank1; 
% change1 = drank1; change2 = drank1; change3 = drank1;
% v = w*Xf';
% [whatever, rank] = sort(v,'descend');

for iter = 1:maxiter
    % get w estimation from standard conjoint
    Z = zeros(nc*nq,30);
    a = 1:size(Xf,1);
    count = 1;
    while count<nc*nq
        id1 = floor(rand*length(a))+1;
        obj1 = a(id1);
        a(id1) = [];
        id2 = floor(rand*length(a))+1;
        obj2 = a(id2);
        a(id2) = [];
        if w*(Xf(obj1,:)-Xf(obj2,:))'>0
            Z(count,:) = Xf(obj1,:)-Xf(obj2,:);
            count = count + 1;
        elseif w*(Xf(obj1,:)-Xf(obj2,:))'<0
            Z(count,:) = Xf(obj2,:)-Xf(obj1,:);
            count = count + 1;
        end
    end
    % Z(abs(Z)<1e-6)=0;
    % Z = unique(Z,'rows');
    Z = sparse(Z);
    model = train(ones(size(Z,1),1),Z,'-s 0 -c 1e6');
    w1(iter,:) = model.w/norm(model.w)*norm(w);
    cov1 = cov([w',w1(iter,:)']);
    rho1(iter) = cov1(1,2)/sqrt(cov1(1,1))/sqrt(cov1(2,2));

    %%
    % get w estimation from adaptive conjoint
    infoappGGBS.X = Xf;
    infoappGGBS.nx = size(Xf,1);
    infoappGGBS.d = 30;
    infoappGGBS.theta = ones(infoappGGBS.nx,1)/infoappGGBS.nx;
    infoappGGBS.maxnq = nc*nq;
    infoappGGBS.nc = nc;
    infoappGGBS.nq = nq;
    infoappGGBS.candidate = 1:infoappGGBS.nx;
    infoappGGBS.update = false;
    [E_app_set,wh_app,pairs2] = pairwiseLearnApproximate_case(infoappGGBS,[],struct([]),w',0);
    w2(iter,:) = wh_app'/norm(wh_app)*norm(w);
    cov2 = cov([w',w2(iter,:)']);
    rho2(iter) = cov2(1,2)/sqrt(cov2(1,1))/sqrt(cov2(2,2));

    %%
    % get w estimation from object identification using GGBS approximation
    pairs3 = [];
    infoappGGBS.X = Xf;
    infoappGGBS.nx = size(Xf,1);
    infoappGGBS.d = 30;
    infoappGGBS.theta = ones(infoappGGBS.nx,1)/infoappGGBS.nx;
    infoappGGBS.theta(rank(1:1000)) = infoappGGBS.theta(rank(1:1000)) + 100*exp(-0.1*(0:999)');
    infoappGGBS.norm = sum(infoappGGBS.theta);
    infoappGGBS.theta = infoappGGBS.theta/sum(infoappGGBS.theta);
    infoappGGBS.tau = 100;
%     infoappGGBS.norm = 1;
    infoappGGBS.maxnq = nc*nq;
    infoappGGBS.nc = nc;
    infoappGGBS.nq = nq;
    infoappGGBS.candidate = 1:infoappGGBS.nx;
    infoappGGBS.update = true;
    [E_app_set,wh_app,pairs3] = pairwiseLearnApproximate_case(infoappGGBS,[],struct([]),w',0);

    Z = zeros(size(pairs3,1),infoappGGBS.d);
    for i = 1:size(pairs3,1)
        Z(i,:) = infoappGGBS.X(pairs3(i,1),:)-infoappGGBS.X(pairs3(i,2),:);
    end
    Z = sparse(Z);
    model = train(ones(size(Z,1),1),Z,'-s 0 -c 1e6');
    w3(iter,:) = model.w/norm(model.w)*norm(w);
    cov3 = cov([w',w3(iter,:)']);
    rho3(iter) = cov3(1,2)/sqrt(cov3(1,1))/sqrt(cov3(2,2));

    %%    
    v = w*Xf';
    [whatever, rank] = sort(v,'descend');
    v1 = w1(iter,:)*Xf';
    [whatever, rank1] = sort(v1,'descend');
    v2 = w2(iter,:)*Xf';
    [whatever, rank2] = sort(v2,'descend');
    v3 = w3(iter,:)*Xf';
    [whatever, rank3] = sort(v3,'descend');
    for i = 1:10
        drank1(iter) = drank1(iter) + abs(find(rank1==rank(i))-i);
        drank2(iter) = drank2(iter) + abs(find(rank2==rank(i))-i);
        drank3(iter) = drank3(iter) + abs(find(rank3==rank(i))-i);
    end

%% design optimization
% design parameters
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
A = zeros(6,15);
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
lb = [0.125*ones(5,1);1;0.5;1;0.25;0.5;0.5;1;1;1;10];
ub = [36;36;24;24;36;200;12;36;24;1.9;1.9;36;36;36;30];
x0 = [2.3;8.87;1.34;1.75;0.41;95.7;0.5;7.44;0.25;0.5;1.9;9.34;11.54;11.57;26.41];

% use newton
% [x,fval,exitflag] = fmincon(@(x)scaleObj(x,y,w),x0,A,b,[],[],lb,ub,@(x)scaleCon(x,y));
% [x1,fval1,exitflag1] = fmincon(@(x)scaleObj(x,y,w1),x0,A,b,[],[],lb,ub,@(x)scaleCon(x,y));
% [x2,fval2,exitflag2] = fmincon(@(x)scaleObj(x,y,w2),x0,A,b,[],[],lb,ub,@(x)scaleCon(x,y));
% [x3,fval3,exitflag3] = fmincon(@(x)scaleObj(x,y,w3),x0,A,b,[],[],lb,ub,@(x)scaleCon(x,y));

% use GA
% options = gaoptimset('InitialPopulation',x);
% [x,fval,exitflag] = ga(@(x)scaleObj(x',y,w),15,A,b,[],[],lb,ub,@(x)scaleCon(x',y),options);
% load trueSolution.mat;
% options = gaoptimset('InitialPopulation',x);
% fval1=0;
% fval2=0;
% fval3=0;
% for i = 1:1
%     [tx,tfval,exitflag] = ga(@(x)scaleObj(x',y,w1(iter,:)),15,A,b,[],[],lb,ub,@(x)scaleCon(x',y),options);
%     if fval1>tfval && exitflag==1
%         x1 = tx;
%         fval1 = tfval;
%         options = gaoptimset('InitialPopulation',x1);
%     end
%     [tx,tfval,exitflag] = ga(@(x)scaleObj(x',y,w2(iter,:)),15,A,b,[],[],lb,ub,@(x)scaleCon(x',y),options);
%     if fval2>tfval && exitflag==1
%         x2 = tx;
%         fval2 = tfval;
%         options = gaoptimset('InitialPopulation',x2);
%     end
%     [tx,tfval,exitflag] = ga(@(x)scaleObj(x',y,w3(iter,:)),15,A,b,[],[],lb,ub,@(x)scaleCon(x',y),options);
%     if fval3>tfval && exitflag==1
%         x3 = tx;
%         fval3 = tfval;
%         options = gaoptimset('InitialPopulation',x3);
%     end
% end
% change1(iter) = fval - scaleObj(x1',y,w);
% change2(iter) = fval - scaleObj(x2',y,w);
% change3(iter) = fval - scaleObj(x3',y,w);

load marketSolution.mat;
% options = gaoptimset('InitialPopulation',x);
% [x,fval,exitflag] = ga(@(x)-calMarketShare(x',y,w),15,A,b,[],[],lb,ub,@(x)scaleCon(x',y),options);
options = gaoptimset('InitialPopulation',x);
fval1=0;
fval2=0;
fval3=0;
for i = 1:5
    [tx,tfval,exitflag] = ga(@(x)-calMarketShare(x',y,w1(iter,:)),15,A,b,[],[],lb,ub,@(x)scaleCon(x',y),options);
    if fval1>tfval && exitflag==1
        x1 = tx;
        fval1 = tfval;
        options = gaoptimset('InitialPopulation',x1);
    end
    [tx,tfval,exitflag] = ga(@(x)-calMarketShare(x',y,w2(iter,:)),15,A,b,[],[],lb,ub,@(x)scaleCon(x',y),options);
    if fval2>tfval && exitflag==1
        x2 = tx;
        fval2 = tfval;
        options = gaoptimset('InitialPopulation',x2);
    end
    [tx,tfval,exitflag] = ga(@(x)-calMarketShare(x',y,w3(iter,:)),15,A,b,[],[],lb,ub,@(x)scaleCon(x',y),options);
    if fval3>tfval && exitflag==1
        x3 = tx;
        fval3 = tfval;
        options = gaoptimset('InitialPopulation',x3);
    end
end
change1(iter) = calMarketShare(x',y,w) - calMarketShare(x1',y,w);
change2(iter) = calMarketShare(x',y,w) - calMarketShare(x2',y,w);
change3(iter) = calMarketShare(x',y,w) - calMarketShare(x3',y,w);

end
allchange = abs([change1,change2,change3]);
[h,pvalue]=ttest2(allchange(:,1),allchange(:,3),0.05,'right','unequal')
% save(saveto);