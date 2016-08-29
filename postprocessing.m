addpath('C:\Users\Max Yi Ren\Documents\Research\Code\Tools\');
addpath('C:\Users\Max Yi Ren\Documents\Research\Code\Tools\liblinear\matlab');
NQ = 50;
TEST = 12;
load basedata.mat;

%% check w convergence
% calculate w under all information
% A = bsxfun(@times,dX,2*(dX*w'>0)-1);
% A(abs(dX*w')<1e-3,:) = [];
% addpath('C:\Users\Max Yi Ren\Documents\Research\Code\Tools\liblinear\matlab');
% model = train([ones(size(A,1),1);-ones(size(A,1),1)],...
%                     sparse([A;-A]),'-s 0 -c 1e6 -e 1e-6 -q');
% w_all = model.w;
load w_all.mat;
DX = dXID + dXID';

% calculate w under ggbs approach
load profit_s100000_inq10.mat
dist_w_ggbs = zeros(TEST,NQ);

for i = 1:TEST
    pairi = pairs_set{i};
    for j = 1:NQ
        Ai = dX(DX(sub2ind(size(DX),pairi(1:j,1),pairi(1:j,2))),:);
        Ai = bsxfun(@times,Ai,sign(pairi(1:j,2)-pairi(1:j,1)));
        modeli = train([ones(size(Ai,1),1);-ones(size(Ai,1),1)],...
                        sparse([Ai;-Ai]),'-s 0 -c 1e6 -e 1e-6 -q');
        wi = modeli.w;
        dist_w_ggbs(i,j) = wi/norm(wi)*w'/norm(w);
    end
end

% calculate w under ggbs approach w/ 6 designs
load profit_s100000_inq6.mat
dist_w_ggbs6 = zeros(TEST,NQ);

for i = 1:TEST
    pairi = pairs_set{i};
    for j = 1:NQ
        Ai = dX(DX(sub2ind(size(DX),pairi(1:j,1),pairi(1:j,2))),:);
        Ai = bsxfun(@times,Ai,sign(pairi(1:j,2)-pairi(1:j,1)));
        modeli = train([ones(size(Ai,1),1);-ones(size(Ai,1),1)],...
                        sparse([Ai;-Ai]),'-s 0 -c 1e6 -e 1e-6 -q');
        wi = modeli.w;
        dist_w_ggbs6(i,j) = wi/norm(wi)*w'/norm(w);
    end
end

% calculate w under ggbs approach w/ 8 designs
load profit_s100000_inq8.mat
dist_w_ggbs8 = zeros(TEST,NQ);

for i = 1:TEST
    pairi = pairs_set{i};
    for j = 1:NQ
        Ai = dX(DX(sub2ind(size(DX),pairi(1:j,1),pairi(1:j,2))),:);
        Ai = bsxfun(@times,Ai,sign(pairi(1:j,2)-pairi(1:j,1)));
        modeli = train([ones(size(Ai,1),1);-ones(size(Ai,1),1)],...
                        sparse([Ai;-Ai]),'-s 0 -c 1e6 -e 1e-6 -q');
        wi = modeli.w;
        dist_w_ggbs8(i,j) = wi/norm(wi)*w'/norm(w);
    end
end


% calculate w under GBS approach
load gbs_s100000_inq10.mat
dist_w_gbs = zeros(TEST,NQ);

for i = 1:TEST
    pairi = pairs_set{i};
    for j = 1:NQ
        Ai = dX(DX(sub2ind(size(DX),pairi(1:j,1),pairi(1:j,2))),:);
        Ai = bsxfun(@times,Ai,sign(pairi(1:j,2)-pairi(1:j,1)));
        modeli = train([ones(size(Ai,1),1);-ones(size(Ai,1),1)],...
                        sparse([Ai;-Ai]),'-s 0 -c 1e6 -e 1e-6 -q');
        wi = modeli.w;
        dist_w_gbs(i,j) = wi/norm(wi)*w'/norm(w);
    end
end

% calculate w under bestprob approach
load bestprob_s100000_inq10.mat
dist_w_bestprob = zeros(1,NQ);

for i = 1:1
    pairi = pairs_set{i};
    for j = 1:NQ
        Ai = dX(DX(sub2ind(size(DX),pairi(1:j,1),pairi(1:j,2))),:);
        Ai = bsxfun(@times,Ai,sign(pairi(1:j,2)-pairi(1:j,1)));
        modeli = train([ones(size(Ai,1),1);-ones(size(Ai,1),1)],...
                        sparse([Ai;-Ai]),'-s 0 -c 1e6 -e 1e-6 -q');
        wi = modeli.w;
        dist_w_bestprob(i,j) = wi/norm(wi)*w'/norm(w);
    end
end

% calculate w under Toubia approach
load toubia_s100000_inq10.mat
dist_w_toubia = zeros(1,NQ);

for i = 1:1
    pairi = pairs_set{i};
    for j = 1:NQ
        Ai = dX(DX(sub2ind(size(DX),pairi(1:j,1),pairi(1:j,2))),:);
        Ai = bsxfun(@times,Ai,sign(pairi(1:j,2)-pairi(1:j,1)));
        modeli = train([ones(size(Ai,1),1);-ones(size(Ai,1),1)],...
                        sparse([Ai;-Ai]),'-s 0 -c 1e6 -e 1e-6 -q');
        wi = modeli.w;
        dist_w_toubia(i,j) = wi/norm(wi)*w'/norm(w);
    end
end

%% check JS divergence
target = sparse(zeros(2455,1));
target(1704) = 0.3709;
target(1705) = 0.6291;

% ggbs approach
load profit_s100000_inq10.mat
jsdiv_ggbs = zeros(TEST,NQ);
found_ggbs = zeros(TEST,1);
for i = 1:TEST
    prob = prob_set{i};
    for j = 1:NQ
        jsdiv_ggbs(i,j) = JSdivergence(target,prob(:,j));
        sorted_prob = sort(prob(:,j),'descend');
        if(prob(1704,j)>=sorted_prob(2))&&found_ggbs(i)==0
            found_ggbs(i) = j;
        end
    end
end

% ggbs approach w/ 6 designs
load profit_s100000_inq6.mat
jsdiv_ggbs6 = zeros(TEST,NQ);
found_ggbs6 = zeros(TEST,1);
for i = 1:TEST
    prob = prob_set{i};
    for j = 1:NQ
        jsdiv_ggbs6(i,j) = JSdivergence(target,prob(:,j));
        sorted_prob = sort(prob(:,j),'descend');
        if(prob(1704,j)>=sorted_prob(2))&&found_ggbs6(i)==0
            found_ggbs6(i) = j;
        end
    end
end

% ggbs approach w/ 8 designs
load profit_s100000_inq8.mat
jsdiv_ggbs8 = zeros(TEST,NQ);
found_ggbs8 = zeros(TEST,1);
for i = 1:TEST
    prob = prob_set{i};
    for j = 1:NQ
        jsdiv_ggbs8(i,j) = JSdivergence(target,prob(:,j));
        sorted_prob = sort(prob(:,j),'descend');
        if(prob(1704,j)>=sorted_prob(2))&&found_ggbs8(i)==0
            found_ggbs8(i) = j;
        end
    end
end

% gbs approach
load gbs_s100000_inq10.mat
jsdiv_gbs = zeros(TEST,NQ);
found_gbs = zeros(TEST,1);
for i = 1:TEST
    prob = prob_set{i};
    for j = 1:NQ
        jsdiv_gbs(i,j) = JSdivergence(target,prob(:,j));
        sorted_prob = sort(prob(:,j),'descend');
        if(prob(1704,j)>=sorted_prob(2))&&found_gbs(i)==0
            found_gbs(i) = j;
        end
    end
end

% bestprob approach
load bestprob_s100000_inq10.mat
jsdiv_bestprob = zeros(1,NQ);
found_bestprob = zeros(1,1);
for i = 1:1
    prob = prob_set{i};
    for j = 1:NQ
        jsdiv_bestprob(i,j) = JSdivergence(target,prob(:,j));
        sorted_prob = sort(prob(:,j),'descend');
        if(prob(1704,j)>=sorted_prob(2))&&found_bestprob(i)==0
            found_bestprob(i) = j;
        end
    end
end

% toubia approach
load toubia_s100000_inq10.mat
jsdiv_toubia = zeros(1,NQ);
found_toubia = zeros(1,1);
for i = 1:1
    prob = prob_set{i};
    for j = 1:NQ
        jsdiv_toubia(i,j) = JSdivergence(target,prob(:,j));
        sorted_prob = sort(prob(:,j),'descend');
        if(prob(1704,j)>=sorted_prob(2))&&found_toubia(i)==0
            found_toubia(i) = j;
        end
    end
end
%% plot
figure; hold on;
shadedErrorBar(1:NQ,jsdiv_ggbs,{@mean,@std},'k',1);
shadedErrorBar(1:NQ,jsdiv_ggbs6,{@mean,@std},'k',1);
shadedErrorBar(1:NQ,jsdiv_ggbs8,{@mean,@std},'k',1);

% plot(jsdiv_ggbs(1,:),'k');
% shadedErrorBar(1:NQ,jsdiv_gbs,{@mean,@std},'r',1);
% shadedErrorBar(1:NQ,jsdiv_bestprob,{@mean,@std},'g',1);
% plot(jsdiv_bestprob(1,:),'g');
% shadedErrorBar(1:NQ,jsdiv_toubia,{@mean,@std},'g',1);
plot(jsdiv_toubia(1,:),'y');

figure; hold on;
shadedErrorBar(1:NQ,dist_w_ggbs,{@mean,@std},'k',1);
shadedErrorBar(1:NQ,dist_w_ggbs6,{@mean,@std},'k',1);
shadedErrorBar(1:NQ,dist_w_ggbs8,{@mean,@std},'k',1);
% plot(dist_w_ggbs(1,:),'k');
% shadedErrorBar(1:NQ,dist_w_gbs,{@mean,@std},'r',1);
% shadedErrorBar(1:NQ,dist_w_bestprob,{@mean,@std},'g',1);
% plot(dist_w_bestprob(1,:),'g');
% shadedErrorBar(1:NQ,dist_w_toubia(1,:),{@mean,@std},'g',1);
plot(dist_w_toubia(1,:),'y');