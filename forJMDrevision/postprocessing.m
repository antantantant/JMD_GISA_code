%% head
% load('..\basedata.mat');
% addpath('..\..\Tools\liblinear\matlab');
% DX = dXID + dXID';
TEST = 12;
target_best_set = 1701*ones(TEST,1); % for the current test, all have 1701 as the target

%% Plots for GGBS
load('profit_s10000_inq100_n1000_comp1_0428.mat');
MAX_ITER = 1000;
prob_ggbs = zeros(TEST,MAX_ITER);
best_ggbs = zeros(TEST,MAX_ITER);
corr_ggbs = zeros(TEST,MAX_ITER);
dist_ggbs = zeros(TEST,MAX_ITER);

for i = 1:TEST
    prob = prob_set{i};
    prob_ggbs(i,:) = prob(target_best_set(i),1:MAX_ITER);
    best_ggbs(i,:) = prob(target_best_set(i),1:MAX_ITER)==max(prob(:,1:MAX_ITER));
    corr_ggbs(i,:) = corr(partworths_set{i},wtrue')';
    dist_ggbs(i,:) = sqrt(sum(bsxfun(@minus,partworths_set{i},wtrue').^2,1));
%     
%     pp = partworths_set{i};
%     for j = 2:1001
%         if j == 220
%             wait = 1;
%         end
%         if (2*norm(pp(:,j))<norm(pp(:,j-1)))
%             pp(:,j) =  pp(:,j)/pp(1,j)*pp(1,j-1);
%             dist_ggbs(i,j) = dist_ggbs(i,j-1);
%             corr_ggbs(i,j) = corr_ggbs(i,j-1);
%         end
%     end
end

% figure; hold on;
% plot(prob_ggbs');
% plot(mean(prob_ggbs,1),'r');
% plot(mean(best_ggbs))

%% Plots for Toubia's method
load('toubia_s10000_n1000_comp1_theta1000_0430.mat');
MAX_ITER = 1000;
prob_toubia = zeros(TEST,MAX_ITER);
corr_toubia = zeros(TEST,MAX_ITER);
dist_toubia = zeros(TEST,MAX_ITER);
best_toubia = zeros(TEST,MAX_ITER);
for i = 1:TEST
    prob = prob_set{i};
    prob_toubia(i,:) = prob(target_best_set(i),1:MAX_ITER);
    best_toubia(i,:) = prob(target_best_set(i),1:MAX_ITER)==max(prob(:,1:MAX_ITER));
    partworths = partworths_set{i};
    corr_toubia(i,:) = corr(partworths,wtrue')';
    dist_toubia(i,:) = sqrt(sum(bsxfun(@minus,partworths,wtrue').^2,1));
end

figure;
subplot(4,1,1);
hold on;
% plot(mean(best_ggbs,1),'-r','LineWidth',1);
plot(mean(prob_ggbs,1),'r','LineWidth',2);
% plot(mean(best_toubia,1),'-b','LineWidth',1);
plot(mean(prob_toubia,1),'b','LineWidth',2);
xlim([0 1000]);
set(gca,'FontSize',16,'Fontname','Timesnewroman');
ylhand = get(gca,'ylabel');
set(ylhand,'string','\pi^{(k^*)}_a','fontsize',20,'Fontname','Timesnewroman');

subplot(4,1,2);
hold on;
plot(mean(best_ggbs,1),'r','LineWidth',2);
plot(mean(best_toubia,1),'b','LineWidth',2);
xlim([0 1000]);
set(gca,'fontSize',16,'fontname','timesnewroman');
ylhand = get(gca,'ylabel');
set(ylhand,'string','best','fontsize',16,'fontname','timesnewroman');
xlhand = get(gca,'xlabel');
set(xlhand,'string','number of queries','fontsize',16,'fontname','timesnewroman');

subplot(4,1,3);
hold on;
plot(mean(corr_ggbs,1),'r','LineWidth',2);
plot(mean(corr_toubia,1),'b','LineWidth',2);
xlim([0 1000]);
set(gca,'fontSize',16,'fontname','timesnewroman');
ylhand = get(gca,'ylabel');
set(ylhand,'string','corr({\bf w}^*,{\bf w}_0)','fontsize',16,'fontname','timesnewroman');
xlhand = get(gca,'xlabel');
set(xlhand,'string','number of queries','fontsize',16,'fontname','timesnewroman');

subplot(4,1,4);
hold on;
plot(mean(dist_ggbs,1),'r','LineWidth',2);
plot(mean(dist_toubia,1),'b','LineWidth',2);
ylabel('||{\bf w}^*-{\bf w}_0||');
xlabel('number of queries');

% savefig('C10q100p.fig');


%% long term
load('profit_s100000_inq100_C1_0210.mat');
TEST = length(prob_set);
prob_ggbs = zeros(TEST,10001);
best_ggbs = zeros(TEST,10001);
for i = 1:TEST
    prob = prob_set{i};
    prob_ggbs(i,:) = prob(1704,1:10001);
    best_ggbs(i,:) = prob(1704,1:10001)==max(prob(:,1:10001));
end

load('toubia_s100000_C1_n10000_0220.mat');
prob_toubia = zeros(TEST,10001);
corr_toubia = zeros(TEST,10001);
for i = 1:TEST
    prob = prob_set{i};
%     pairs = pairs_set{i};
    prob_toubia(i,:) = prob(1704,1:10001);
    best_toubia(i,:) = prob(1704,1:10001)==max(prob(:,1:10001));
end
figure;
hold on;
shadedErrorBar(1:10001,prob_ggbs,{@mean,@std},'k',1);
shadedErrorBar(1:10001,prob_toubia,{@mean,@std},'k',1);

plot(mean(best_ggbs,1),'k','LineWidth',1);
% plot(mean(prob_ggbs,1),'k','LineWidth',2);
plot(mean(best_toubia,1),'-.k','LineWidth',1);
% plot(mean(prob_toubia,1),'-.k','LineWidth',2);
ylabel('\pi_{\Theta^a_{k^*}}');

%% change candidate query size
load('profit_s100000_inq100_C1_n1000_0224.mat');
TEST = length(prob_set);
prob_ggbs_q100 = zeros(TEST,1001);
best_ggbs_q100 = zeros(TEST,1001);
for i = 1:TEST
    prob = prob_set{i};
    prob_ggbs_q100(i,:) = prob(1704,1:1001);
    best_ggbs_q100(i,:) = prob(1704,1:1001)==max(prob(:,1:1001));
end

load('profit_s100000_inq1000_C1_n1000_0222.mat');
TEST = length(prob_set);
prob_ggbs_q1000 = zeros(TEST,1001);
best_ggbs_q1000 = zeros(TEST,1001);
for i = 1:TEST
    prob = prob_set{i};
    prob_ggbs_q1000(i,:) = prob(1704,1:1001);
    best_ggbs_q1000(i,:) = prob(1704,1:1001)==max(prob(:,1:1001));
end

load('profit_s100000_inq10_C1_n1000_0222.mat');
TEST = length(prob_set);
prob_ggbs_q10 = zeros(TEST,1001);
best_ggbs_q10 = zeros(TEST,1001);
for i = 1:TEST
    prob = prob_set{i};
    prob_ggbs_q10(i,:) = prob(1704,1:1001);
    best_ggbs_q10(i,:) = prob(1704,1:1001)==max(prob(:,1:1001));
end

figure;
hold on;
shadedErrorBar(1:1001,prob_ggbs_q10,{@mean,@std},'g',1);
shadedErrorBar(1:1001,prob_ggbs_q100,{@mean,@std},'r',1);
shadedErrorBar(1:1001,prob_ggbs_q1000,{@mean,@std},'b',1);
plot(mean(best_ggbs_q10,1),'g','LineWidth',2);
plot(mean(best_ggbs_q100,1),'r','LineWidth',2);
plot(mean(best_ggbs_q1000,1),'b','LineWidth',2);
xlabel('number of queries');
ylabel('\pi_{\Theta^a_{k^*}}');

%% debug
load('profit_s100000_inq100_C1_n1000_0224.mat');
TEST = length(prob_set);
prob_ggbs_0224 = zeros(TEST,1001);
best_ggbs_0224 = zeros(TEST,1001);
for i = 1:TEST
    prob = prob_set{i};
    prob_ggbs_0224(i,:) = prob(1704,1:1001);
    best_ggbs_0224(i,:) = prob(1704,1:1001)==max(prob(:,1:1001));
end

load('profit_s100000_inq100_C1_0210.mat');
TEST = length(prob_set);
prob_ggbs_0210 = zeros(TEST,1001);
best_ggbs_0210 = zeros(TEST,1001);
for i = 1:TEST
    prob = prob_set{i};
    prob_ggbs_0210(i,:) = prob(1704,1:1001);
    best_ggbs_0210(i,:) = prob(1704,1:1001)==max(prob(:,1:1001));
end

figure;
hold on;
shadedErrorBar(1:1001,prob_ggbs_0224,{@mean,@std},'r',1);
shadedErrorBar(1:1001,prob_ggbs_0210,{@mean,@std},'b',1);

%% with and without Abernethy
load('profit_s100000_inq100_C1_n1000_0224.mat');
TEST = length(prob_set);
prob_ggbs_0224 = zeros(TEST,1001);
best_ggbs_0224 = zeros(TEST,1001);
for i = 1:TEST
    prob = prob_set{i};
    prob_ggbs_0224(i,:) = prob(1704,1:1001);
    best_ggbs_0224(i,:) = prob(1704,1:1001)==max(prob(:,1:1001));
end

load('profit_s100000_inq100_C1_n1000_0228.mat');
TEST = length(prob_set);
prob_ggbs_0228 = zeros(TEST,1001);
best_ggbs_0228 = zeros(TEST,1001);
for i = 1:TEST
    prob = prob_set{i};
    prob_ggbs_0228(i,:) = prob(1704,1:1001);
    best_ggbs_0228(i,:) = prob(1704,1:1001)==max(prob(:,1:1001));
end

figure;
hold on;
shadedErrorBar(1:1001,prob_ggbs_0224,{@mean,@std},'r',1);
shadedErrorBar(1:1001,prob_ggbs_0228,{@mean,@std},'b',1);
shadedErrorBar(1:1001,prob_toubia(:,1:1001),{@mean,@std},'g',1);


%% investigate C=0.1 Toubia corr->1 dist->0 but probability is not close to 1
load('toubia_s100000_C01_n1000_0218.mat');
C = 0.1;
pairs = pairs_set{1};
A = dX(DX(sub2ind(size(DX),pairs(:,1),pairs(:,2))),:);
A = bsxfun(@times,A,sign(pairs(:,2)-pairs(:,1)));
% b = zeros(size(pairs,1),1);
% model = train(ones(size(A,1),1),...
%                     sparse(A), sprintf('-s 0 -e 1e-6 -q -c %f', C));
% w0_toubia = model.w;
ww = partworths_set{1};
w0_toubia = ww(:,end)';
As = bsxfun(@times, A, sqrt(exp(A*w0_toubia'))./(1+exp(A*w0_toubia')));
Sigma_inv = (eye(d)/C+(As'*As));
W = W0*(inv(chol(Sigma_inv)))';
W_toubia = bsxfun(@plus, W, w0_toubia);
util = (W_toubia*Xf');
obj_app = bsxfun(@plus,-log(1+exp(-util)),log(c'));
y_toubia = obj_app(:,1704)==max(obj_app,[],2);

load('profit_s100000_inq100_C01_n1000_0226.mat');
C = 0.1;
pairs = pairs_set{1};
A = dX(DX(sub2ind(size(DX),pairs(:,1),pairs(:,2))),:);
A = bsxfun(@times,A,sign(pairs(:,2)-pairs(:,1)));
% b = zeros(size(pairs,1),1);
% model = train(ones(size(A,1),1),...
%                     sparse(A), sprintf('-s 0 -e 1e-6 -q -c %f', C));
% w0_gisa = model.w;
ww = partworths_set{1};
w0_gisa = ww(:,end)';
As = bsxfun(@times, A, sqrt(exp(A*w0_gisa'))./(1+exp(A*w0_gisa')));
Sigma_inv = (eye(d)/C+(As'*As));
W = W0*(inv(chol(Sigma_inv)))';
W_gisa = bsxfun(@plus, W, w0_gisa);
util = (W_gisa*Xf');
obj_app = bsxfun(@plus,-log(1+exp(-util)),log(c'));
y_gisa = obj_app(:,1704)==max(obj_app,[],2);

%project to w0_gisa-w0_toubia
W_gisa_p = W_gisa(y_gisa,:);
W_gisa_n = W_gisa(~y_gisa,:);
v = w0_gisa-w0_toubia;
v = v/norm(v);
pW_gisa_p = W_gisa_p*v';
pW_gisa_n = W_gisa_n*v';

W_toubia_p = W_toubia(y_toubia,:);
W_toubia_n = W_toubia(~y_toubia,:);
pW_toubia_p = W_toubia_p*v';
pW_toubia_n = W_toubia_n*v';

pw = w*v';

figure;hold on;
data.h4 = [pW_toubia_n;pW_toubia_p];
data.h2 = [pW_gisa_n;pW_gisa_p];
data.h1 = pW_gisa_p;
data.h3 = pW_toubia_p;
% h4 = nhist([pW_toubia_n;pW_toubia_p]);
% h2 = nhist([pW_gisa_n;pW_gisa_p]);
% h1 = nhist(pW_gisa_p);
% h3 = nhist(pW_toubia_p);
nhist(data,'numbers');


%random projection but perpendicular to ...
W_gisa_p = W_gisa(y_gisa,:);
W_gisa_n = W_gisa(~y_gisa,:);
v = rand(30,1);
dw = mean(W_gisa_p)-mean(W_gisa_n);
v(30) = -(dw(1:29)*v(1:29))/dw(30);
v = v/norm(v);
v = [v,dw'/norm(dw)];
pw = w*v;
pW_gisa = W_gisa*v;
figure; hold on;
% plot(pW_gisa(y_gisa,1),pW_gisa(y_gisa,2),'sr','MarkerSize',4);
% plot(pW_gisa(~y_gisa,1),pW_gisa(~y_gisa,2),'.b','MarkerSize',8);
% plot(pw(1),pw(2),'.k','MarkerSize',20);
plot(W_gisa(y_gisa,1),W_gisa(y_gisa,2),'sr','MarkerSize',4);
plot(W_gisa(~y_gisa,1),W_gisa(~y_gisa,2),'.b','MarkerSize',8);
plot(pw(1),pw(2),'.k','MarkerSize',20);

W_toubia_p = W_toubia(y_toubia,:);
W_toubia_n = W_gisa(~y_toubia,:);
v = rand(30,1);
dw = mean(W_toubia_p)-mean(W_toubia_n);
v(30) = -(dw(1:29)*v(1:29))/dw(30);
v = v/norm(v);
v = [v,dw'/norm(dw)];
pw = w*v;
pW_toubia = W_toubia*v;
figure; hold on;
plot(pW_toubia(y_toubia,1),pW_gisa(y_toubia,2),'sr','MarkerSize',4);
plot(pW_toubia(~y_toubia,1),pW_gisa(~y_toubia,2),'.b','MarkerSize',8);
plot(pw(1),pw(2),'.k','MarkerSize',20);
