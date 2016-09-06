%% Produce figures for the journal paper
%% head
load('../basedata.mat');
% addpath('..\..\Tools\liblinear\matlab');
DX = dXID + dXID';
TEST = 20;
MAX_ITER = 1000;
inq = 100;
theta = 100;
s = 1e4;

c = Xf(:,26:30)*price'-cv; %price - cost

XID = 1:30; 
XID(5:5:30)=[];

nt = size(Xf,1); % number of testing object

wtrue = w*theta;
wtrue(1:5) = wtrue(1:5)-wtrue(5);
wtrue(6:10) = wtrue(6:10)-wtrue(10);
wtrue(11:15) = wtrue(11:15)-wtrue(15);
wtrue(16:20) = wtrue(16:20)-wtrue(20);
wtrue(21:25) = wtrue(21:25)-wtrue(25);
wtrue(26:30) = wtrue(26:30)-wtrue(30);
num_competitor = 1;

%% Plots for GISA and Abernethy's methods
plotGISA = true;
plotAbernethy = true;
theta = 1;
s = 1e3;
if(plotGISA)
    load(['gisa_s',num2str(s),'_inq',num2str(inq),'_n',num2str(MAX_ITER),...
    '_comp',num2str(num_competitor),'_theta',num2str(theta),...
    '_nt',num2str(nt),'_09042016.mat']);
    prob_gisa = zeros(TEST,MAX_ITER);
    best_gisa = zeros(TEST,MAX_ITER);
    corr_gisa = zeros(TEST,MAX_ITER);
    dist_gisa = zeros(TEST,MAX_ITER);
    prob_gisa_b = zeros(100,MAX_ITER);
    best_gisa_b = zeros(100,MAX_ITER);
    corr_gisa_b = zeros(100,MAX_ITER);
    dist_gisa_b = zeros(100,MAX_ITER);
    for i = 1:TEST
        prob = prob_set{i};
        prob_gisa(i,:) = prob(target_best_set(i),1:MAX_ITER);
        best_gisa(i,:) = prob(target_best_set(i),1:MAX_ITER)==max(prob(:,1:MAX_ITER));
        corr_gisa(i,:) = corr(partworths_set{i},wtrue')';
        dist_gisa(i,:) = sqrt(sum(bsxfun(@minus,partworths_set{i},wtrue').^2,1));
    end
    for i = 1:MAX_ITER
        prob_gisa_b(:,i) = bootstrp(100,@mean,prob_gisa(:,i));
        best_gisa_b(:,i) = bootstrp(100,@mean,best_gisa(:,i));
        corr_gisa_b(:,i) = bootstrp(100,@mean,corr_gisa(:,i));
        dist_gisa_b(:,i) = bootstrp(100,@mean,dist_gisa(:,i));
    end
end
s = 1e4;
if(plotAbernethy)
    load(['abernethy_s',num2str(s),'_n',num2str(MAX_ITER),...
    '_comp',num2str(num_competitor),'_theta',num2str(theta),...
    '_nt',num2str(nt),'_09042016.mat']);
    prob_abernethy = zeros(TEST,MAX_ITER);
    corr_abernethy = zeros(TEST,MAX_ITER);
    dist_abernethy = zeros(TEST,MAX_ITER);
    best_abernethy = zeros(TEST,MAX_ITER);
    prob_abernethy_b = zeros(100,MAX_ITER);
    corr_abernethy_b = zeros(100,MAX_ITER);
    dist_abernethy_b = zeros(100,MAX_ITER);
    best_abernethy_b = zeros(100,MAX_ITER);
    for i = 1:TEST
        prob = prob_set{i};
        prob_abernethy(i,:) = prob(target_best_set(i),1:MAX_ITER);
        best_abernethy(i,:) = prob(target_best_set(i),1:MAX_ITER)==max(prob(:,1:MAX_ITER));
        partworths = partworths_set{i};
        corr_abernethy(i,:) = corr(partworths,wtrue')';
        dist_abernethy(i,:) = sqrt(sum(bsxfun(@minus,partworths,wtrue').^2,1));
    end
    for i = 1:MAX_ITER
        prob_abernethy_b(:,i) = bootstrp(100,@mean,prob_abernethy(:,i));
        best_abernethy_b(:,i) = bootstrp(100,@mean,best_abernethy(:,i));
        corr_abernethy_b(:,i) = bootstrp(100,@mean,corr_abernethy(:,i));
        dist_abernethy_b(:,i) = bootstrp(100,@mean,dist_abernethy(:,i));
    end
end

figure;

subplot(4,1,1);
hold on;
title(['Abernethy vs.'*plotAbernethy,'GISA, '*plotGISA,'theta=',num2str(theta)]);

if(plotGISA)
    shadedErrorBar(1:MAX_ITER,prob_gisa_b,{@mean,@std},'r',2);
end
if(plotAbernethy)
    shadedErrorBar(1:MAX_ITER,prob_abernethy_b,{@mean,@std},'b',2);
end
% plot(mean(prob_ggbs,1),'r','LineWidth',2);
% plot(mean(prob_toubia,1),'b','LineWidth',2);
xlim([1 MAX_ITER]);
set(gca,'FontSize',16,'Fontname','Timesnewroman');
ylhand = get(gca,'ylabel');
set(ylhand,'string','\pi^{(k^*)}_a','fontsize',20,'Fontname','Timesnewroman');
xlhand = get(gca,'xlabel');
set(xlhand,'string','number of queries','fontsize',16,'fontname','timesnewroman');
plot([1,MAX_ITER],[1,1],'.-k','LineWidth',2);
plot([1,MAX_ITER],[0,0],'.-k','LineWidth',2);

subplot(4,1,2);
hold on;
if(plotGISA)
    shadedErrorBar(1:MAX_ITER,best_gisa_b,{@mean,@std},'r',2);
end
if(plotAbernethy)
    shadedErrorBar(1:MAX_ITER,best_abernethy_b,{@mean,@std},'b',2);
end
% plot(mean(best_ggbs,1),'r','LineWidth',2);
% plot(mean(best_toubia,1),'b','LineWidth',2);
xlim([1 MAX_ITER]);
set(gca,'fontSize',16,'fontname','timesnewroman');
ylhand = get(gca,'ylabel');
set(ylhand,'string','best','fontsize',16,'fontname','timesnewroman');
xlhand = get(gca,'xlabel');
set(xlhand,'string','number of queries','fontsize',16,'fontname','timesnewroman');
plot([1,MAX_ITER],[1,1],'.-k','LineWidth',2);
plot([1,MAX_ITER],[0,0],'.-k','LineWidth',2);

subplot(4,1,3);
hold on;
if(plotGISA)
    shadedErrorBar(1:MAX_ITER,corr_gisa_b,{@mean,@std},'r',2);
end
if(plotAbernethy)
    shadedErrorBar(1:MAX_ITER,corr_abernethy_b,{@mean,@std},'b',2);
end
% plot(mean(corr_ggbs,1),'r','LineWidth',2);
% plot(mean(corr_toubia,1),'b','LineWidth',2);
xlim([1 MAX_ITER]);
set(gca,'fontSize',16,'fontname','timesnewroman');
ylhand = get(gca,'ylabel');
set(ylhand,'string','corr({\bf w}^*,{\bf w}_0)','fontsize',16,'fontname','timesnewroman');
xlhand = get(gca,'xlabel');
set(xlhand,'string','number of queries','fontsize',16,'fontname','timesnewroman');
plot([1,MAX_ITER],[1,1],'.-k','LineWidth',2);
plot([1,MAX_ITER],[0,0],'.-k','LineWidth',2);

subplot(4,1,4);
hold on;
if(plotGISA)
    shadedErrorBar(1:MAX_ITER,dist_gisa_b,{@mean,@std},'r',2);
end
if(plotAbernethy)
    shadedErrorBar(1:MAX_ITER,dist_abernethy_b,{@mean,@std},'b',2);
end
% plot(mean(dist_ggbs,1),'r','LineWidth',2);
% plot(mean(dist_toubia,1),'b','LineWidth',2);
xlim([1 MAX_ITER]);
set(gca,'fontSize',16,'fontname','timesnewroman');
ylhand = get(gca,'ylabel');
set(ylhand,'string','||{\bf w}^*-{\bf w}_0||','fontsize',16,'fontname','timesnewroman');
xlhand = get(gca,'xlabel');
set(xlhand,'string','number of queries','fontsize',16,'fontname','timesnewroman');

% savefig('C10q100p.fig');

%% Plots for GISA s=1000 vs s=10000
plotGISA = true;
plotAbernethy = true;
inq = 100;
s = 1e4;
if(plotGISA)
    load(['gisa_s',num2str(s),'_inq',num2str(inq),'_n',num2str(MAX_ITER),...
    '_comp',num2str(num_competitor),'_theta',num2str(theta),...
    '_nt',num2str(nt),'_09042016.mat']);
    prob_gisa = zeros(TEST,MAX_ITER);
    best_gisa = zeros(TEST,MAX_ITER);
    corr_gisa = zeros(TEST,MAX_ITER);
    dist_gisa = zeros(TEST,MAX_ITER);
    prob_gisa_b = zeros(100,MAX_ITER);
    best_gisa_b = zeros(100,MAX_ITER);
    corr_gisa_b = zeros(100,MAX_ITER);
    dist_gisa_b = zeros(100,MAX_ITER);
    for i = 1:TEST
        prob = prob_set{i};
        prob_gisa(i,:) = prob(target_best_set(i),1:MAX_ITER);
        best_gisa(i,:) = prob(target_best_set(i),1:MAX_ITER)==max(prob(:,1:MAX_ITER));
        corr_gisa(i,:) = corr(partworths_set{i},wtrue')';
        dist_gisa(i,:) = sqrt(sum(bsxfun(@minus,partworths_set{i},wtrue').^2,1));
    end
    for i = 1:MAX_ITER
        prob_gisa_b(:,i) = bootstrp(100,@mean,prob_gisa(:,i));
        best_gisa_b(:,i) = bootstrp(100,@mean,best_gisa(:,i));
        corr_gisa_b(:,i) = bootstrp(100,@mean,corr_gisa(:,i));
        dist_gisa_b(:,i) = bootstrp(100,@mean,dist_gisa(:,i));
    end
end

s = 1e3;
if(plotAbernethy)
    load(['gisa_s',num2str(s),'_inq',num2str(inq),'_n',num2str(MAX_ITER),...
    '_comp',num2str(num_competitor),'_theta',num2str(theta),...
    '_nt',num2str(nt),'_09042016.mat']);
    prob_abernethy = zeros(TEST,MAX_ITER);
    corr_abernethy = zeros(TEST,MAX_ITER);
    dist_abernethy = zeros(TEST,MAX_ITER);
    best_abernethy = zeros(TEST,MAX_ITER);
    prob_abernethy_b = zeros(100,MAX_ITER);
    corr_abernethy_b = zeros(100,MAX_ITER);
    dist_abernethy_b = zeros(100,MAX_ITER);
    best_abernethy_b = zeros(100,MAX_ITER);
    for i = 1:TEST
        prob = prob_set{i};
        prob_abernethy(i,:) = prob(target_best_set(i),1:MAX_ITER);
        best_abernethy(i,:) = prob(target_best_set(i),1:MAX_ITER)==max(prob(:,1:MAX_ITER));
        partworths = partworths_set{i};
        corr_abernethy(i,:) = corr(partworths,wtrue')';
        dist_abernethy(i,:) = sqrt(sum(bsxfun(@minus,partworths,wtrue').^2,1));
    end
    for i = 1:MAX_ITER
        prob_abernethy_b(:,i) = bootstrp(100,@mean,prob_abernethy(:,i));
        best_abernethy_b(:,i) = bootstrp(100,@mean,best_abernethy(:,i));
        corr_abernethy_b(:,i) = bootstrp(100,@mean,corr_abernethy(:,i));
        dist_abernethy_b(:,i) = bootstrp(100,@mean,dist_abernethy(:,i));
    end
end

figure;

subplot(4,1,1);
hold on;
title(['Abernethy vs.'*plotAbernethy,'GISA, '*plotGISA,'theta=',num2str(theta)]);

if(plotGISA)
    shadedErrorBar(1:MAX_ITER,prob_gisa_b,{@mean,@std},'r',2);
end
if(plotAbernethy)
    shadedErrorBar(1:MAX_ITER,prob_abernethy_b,{@mean,@std},'b',2);
end
% plot(mean(prob_ggbs,1),'r','LineWidth',2);
% plot(mean(prob_toubia,1),'b','LineWidth',2);
xlim([1 MAX_ITER]);
set(gca,'FontSize',16,'Fontname','Timesnewroman');
ylhand = get(gca,'ylabel');
set(ylhand,'string','\pi^{(k^*)}_a','fontsize',20,'Fontname','Timesnewroman');
xlhand = get(gca,'xlabel');
set(xlhand,'string','number of queries','fontsize',16,'fontname','timesnewroman');
plot([1,MAX_ITER],[1,1],'.-k','LineWidth',2);
plot([1,MAX_ITER],[0,0],'.-k','LineWidth',2);

subplot(4,1,2);
hold on;
if(plotGISA)
    shadedErrorBar(1:MAX_ITER,best_gisa_b,{@mean,@std},'r',2);
end
if(plotAbernethy)
    shadedErrorBar(1:MAX_ITER,best_abernethy_b,{@mean,@std},'b',2);
end
% plot(mean(best_ggbs,1),'r','LineWidth',2);
% plot(mean(best_toubia,1),'b','LineWidth',2);
xlim([1 MAX_ITER]);
set(gca,'fontSize',16,'fontname','timesnewroman');
ylhand = get(gca,'ylabel');
set(ylhand,'string','best','fontsize',16,'fontname','timesnewroman');
xlhand = get(gca,'xlabel');
set(xlhand,'string','number of queries','fontsize',16,'fontname','timesnewroman');
plot([1,MAX_ITER],[1,1],'.-k','LineWidth',2);
plot([1,MAX_ITER],[0,0],'.-k','LineWidth',2);

subplot(4,1,3);
hold on;
if(plotGISA)
    shadedErrorBar(1:MAX_ITER,corr_gisa_b,{@mean,@std},'r',2);
end
if(plotAbernethy)
    shadedErrorBar(1:MAX_ITER,corr_abernethy_b,{@mean,@std},'b',2);
end
% plot(mean(corr_ggbs,1),'r','LineWidth',2);
% plot(mean(corr_toubia,1),'b','LineWidth',2);
xlim([1 MAX_ITER]);
set(gca,'fontSize',16,'fontname','timesnewroman');
ylhand = get(gca,'ylabel');
set(ylhand,'string','corr({\bf w}^*,{\bf w}_0)','fontsize',16,'fontname','timesnewroman');
xlhand = get(gca,'xlabel');
set(xlhand,'string','number of queries','fontsize',16,'fontname','timesnewroman');
plot([1,MAX_ITER],[1,1],'.-k','LineWidth',2);
plot([1,MAX_ITER],[0,0],'.-k','LineWidth',2);

subplot(4,1,4);
hold on;
if(plotGISA)
    shadedErrorBar(1:MAX_ITER,dist_gisa_b,{@mean,@std},'r',2);
end
if(plotAbernethy)
    shadedErrorBar(1:MAX_ITER,dist_abernethy_b,{@mean,@std},'b',2);
end
% plot(mean(dist_ggbs,1),'r','LineWidth',2);
% plot(mean(dist_toubia,1),'b','LineWidth',2);
xlim([1 MAX_ITER]);
set(gca,'fontSize',16,'fontname','timesnewroman');
ylhand = get(gca,'ylabel');
set(ylhand,'string','||{\bf w}^*-{\bf w}_0||','fontsize',16,'fontname','timesnewroman');
xlhand = get(gca,'xlabel');
set(xlhand,'string','number of queries','fontsize',16,'fontname','timesnewroman');




%% Plots for GISA inq=100 vs inq=10
plotGISA = true;
plotAbernethy = true;
s = 1e4;
inq = 100;
if(plotGISA)
    load(['gisa_s',num2str(s),'_inq',num2str(inq),'_n',num2str(MAX_ITER),...
    '_comp',num2str(num_competitor),'_theta',num2str(theta),...
    '_nt',num2str(nt),'_09042016.mat']);
    prob_gisa = zeros(TEST,MAX_ITER);
    best_gisa = zeros(TEST,MAX_ITER);
    corr_gisa = zeros(TEST,MAX_ITER);
    dist_gisa = zeros(TEST,MAX_ITER);
    prob_gisa_b = zeros(100,MAX_ITER);
    best_gisa_b = zeros(100,MAX_ITER);
    corr_gisa_b = zeros(100,MAX_ITER);
    dist_gisa_b = zeros(100,MAX_ITER);
    for i = 1:TEST
        prob = prob_set{i};
        prob_gisa(i,:) = prob(target_best_set(i),1:MAX_ITER);
        best_gisa(i,:) = prob(target_best_set(i),1:MAX_ITER)==max(prob(:,1:MAX_ITER));
        corr_gisa(i,:) = corr(partworths_set{i},wtrue')';
        dist_gisa(i,:) = sqrt(sum(bsxfun(@minus,partworths_set{i},wtrue').^2,1));
    end
    for i = 1:MAX_ITER
        prob_gisa_b(:,i) = bootstrp(100,@mean,prob_gisa(:,i));
        best_gisa_b(:,i) = bootstrp(100,@mean,best_gisa(:,i));
        corr_gisa_b(:,i) = bootstrp(100,@mean,corr_gisa(:,i));
        dist_gisa_b(:,i) = bootstrp(100,@mean,dist_gisa(:,i));
    end
end

s = 1e3;
inq = 10;
if(plotAbernethy)
    load(['gisa_s',num2str(s),'_inq',num2str(inq),'_n',num2str(MAX_ITER),...
    '_comp',num2str(num_competitor),'_theta',num2str(theta),...
    '_nt',num2str(nt),'_09042016.mat']);
    prob_abernethy = zeros(TEST,MAX_ITER);
    corr_abernethy = zeros(TEST,MAX_ITER);
    dist_abernethy = zeros(TEST,MAX_ITER);
    best_abernethy = zeros(TEST,MAX_ITER);
    prob_abernethy_b = zeros(100,MAX_ITER);
    corr_abernethy_b = zeros(100,MAX_ITER);
    dist_abernethy_b = zeros(100,MAX_ITER);
    best_abernethy_b = zeros(100,MAX_ITER);
    for i = 1:TEST
        prob = prob_set{i};
        prob_abernethy(i,:) = prob(target_best_set(i),1:MAX_ITER);
        best_abernethy(i,:) = prob(target_best_set(i),1:MAX_ITER)==max(prob(:,1:MAX_ITER));
        partworths = partworths_set{i};
        corr_abernethy(i,:) = corr(partworths,wtrue')';
        dist_abernethy(i,:) = sqrt(sum(bsxfun(@minus,partworths,wtrue').^2,1));
    end
    for i = 1:MAX_ITER
        prob_abernethy_b(:,i) = bootstrp(100,@mean,prob_abernethy(:,i));
        best_abernethy_b(:,i) = bootstrp(100,@mean,best_abernethy(:,i));
        corr_abernethy_b(:,i) = bootstrp(100,@mean,corr_abernethy(:,i));
        dist_abernethy_b(:,i) = bootstrp(100,@mean,dist_abernethy(:,i));
    end
end

figure;

subplot(4,1,1);
hold on;
title(['Abernethy vs.'*plotAbernethy,'GISA, '*plotGISA,'theta=',num2str(theta)]);

if(plotGISA)
    shadedErrorBar(1:MAX_ITER,prob_gisa_b,{@mean,@std},'r',2);
end
if(plotAbernethy)
    shadedErrorBar(1:MAX_ITER,prob_abernethy_b,{@mean,@std},'b',2);
end
% plot(mean(prob_ggbs,1),'r','LineWidth',2);
% plot(mean(prob_toubia,1),'b','LineWidth',2);
xlim([1 MAX_ITER]);
set(gca,'FontSize',16,'Fontname','Timesnewroman');
ylhand = get(gca,'ylabel');
set(ylhand,'string','\pi^{(k^*)}_a','fontsize',20,'Fontname','Timesnewroman');
xlhand = get(gca,'xlabel');
set(xlhand,'string','number of queries','fontsize',16,'fontname','timesnewroman');
plot([1,MAX_ITER],[1,1],'.-k','LineWidth',2);
plot([1,MAX_ITER],[0,0],'.-k','LineWidth',2);

subplot(4,1,2);
hold on;
if(plotGISA)
    shadedErrorBar(1:MAX_ITER,best_gisa_b,{@mean,@std},'r',2);
end
if(plotAbernethy)
    shadedErrorBar(1:MAX_ITER,best_abernethy_b,{@mean,@std},'b',2);
end
% plot(mean(best_ggbs,1),'r','LineWidth',2);
% plot(mean(best_toubia,1),'b','LineWidth',2);
xlim([1 MAX_ITER]);
set(gca,'fontSize',16,'fontname','timesnewroman');
ylhand = get(gca,'ylabel');
set(ylhand,'string','best','fontsize',16,'fontname','timesnewroman');
xlhand = get(gca,'xlabel');
set(xlhand,'string','number of queries','fontsize',16,'fontname','timesnewroman');
plot([1,MAX_ITER],[1,1],'.-k','LineWidth',2);
plot([1,MAX_ITER],[0,0],'.-k','LineWidth',2);

subplot(4,1,3);
hold on;
if(plotGISA)
    shadedErrorBar(1:MAX_ITER,corr_gisa_b,{@mean,@std},'r',2);
end
if(plotAbernethy)
    shadedErrorBar(1:MAX_ITER,corr_abernethy_b,{@mean,@std},'b',2);
end
% plot(mean(corr_ggbs,1),'r','LineWidth',2);
% plot(mean(corr_toubia,1),'b','LineWidth',2);
xlim([1 MAX_ITER]);
set(gca,'fontSize',16,'fontname','timesnewroman');
ylhand = get(gca,'ylabel');
set(ylhand,'string','corr({\bf w}^*,{\bf w}_0)','fontsize',16,'fontname','timesnewroman');
xlhand = get(gca,'xlabel');
set(xlhand,'string','number of queries','fontsize',16,'fontname','timesnewroman');
plot([1,MAX_ITER],[1,1],'.-k','LineWidth',2);
plot([1,MAX_ITER],[0,0],'.-k','LineWidth',2);

subplot(4,1,4);
hold on;
if(plotGISA)
    shadedErrorBar(1:MAX_ITER,dist_gisa_b,{@mean,@std},'r',2);
end
if(plotAbernethy)
    shadedErrorBar(1:MAX_ITER,dist_abernethy_b,{@mean,@std},'b',2);
end
% plot(mean(dist_ggbs,1),'r','LineWidth',2);
% plot(mean(dist_toubia,1),'b','LineWidth',2);
xlim([1 MAX_ITER]);
set(gca,'fontSize',16,'fontname','timesnewroman');
ylhand = get(gca,'ylabel');
set(ylhand,'string','||{\bf w}^*-{\bf w}_0||','fontsize',16,'fontname','timesnewroman');
xlhand = get(gca,'xlabel');
set(xlhand,'string','number of queries','fontsize',16,'fontname','timesnewroman');


%% GISA most likely to be the best vs. best expectation
s = 1e4;
inq = 100;

load(['gisa_s',num2str(s),'_inq',num2str(inq),'_n',num2str(MAX_ITER),...
'_comp',num2str(num_competitor),'_theta',num2str(theta),...
'_nt',num2str(nt),'_09042016.mat']);
% prob_gisa = zeros(TEST,MAX_ITER);
% expected_gisa = zeros(TEST,MAX_ITER);
rank_prob_gisa = zeros(TEST,MAX_ITER);
rank_expected_gisa = zeros(TEST,MAX_ITER);

for i = 1:TEST
    prob = prob_set{i};
    expected = expected_value_set{i};
%     prob_gisa(i,:) = prob(target_best_set(i),1:MAX_ITER);
%     expected_gisa(i,:) = expected(target_best_set(i),1:MAX_ITER);
    rank_prob_gisa(i,:) = sum(bsxfun(@ge, prob(:,1:MAX_ITER), prob(target_best_set(i),1:MAX_ITER)),1);
    rank_expected_gisa(i,:) = sum(bsxfun(@ge, expected(:,1:MAX_ITER), expected(target_best_set(i),1:MAX_ITER)),1);
end
figure;
hold on;
title('Rank comparison: Probability to be the best vs. expected objective');

shadedErrorBar(1:MAX_ITER,rank_prob_gisa,{@mean,@std},'r',2);
shadedErrorBar(1:MAX_ITER,rank_expected_gisa,{@mean,@std},'b',2);
xlim([1 MAX_ITER]);
set(gca,'FontSize',16,'Fontname','Timesnewroman');
ylhand = get(gca,'ylabel');
set(ylhand,'string','rank','fontsize',20,'Fontname','Timesnewroman');
xlhand = get(gca,'xlabel');
set(xlhand,'string','number of queries','fontsize',16,'fontname','timesnewroman');
plot([1,MAX_ITER],[nt,nt],'.-k','LineWidth',2);
plot([1,MAX_ITER],[1,1],'.-k','LineWidth',2);
