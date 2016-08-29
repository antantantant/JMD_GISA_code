addpath('C:\Users\Max Yi Ren\Documents\Research\Code\Tools\');
addpath('C:\Users\Max Yi Ren\Documents\Research\Code\Tools\liblinear\matlab');
NQ = 30;
% clear;
% 
% load bestprob_batch0_s1000_inq100_forprofit.mat
% bestprob_t = zeros(20,1);
% for i = 1:20
%     p = prob_set{i};
%     for j = 1:size(p,2)
%         [~,id]=sort(p(:,j),'descend');
%         if id(1)==1705&&id(2)==1704
%             bestprob_t(i) =  j;
%             break;
%         end
%     end
% end
% 
% load ggbs_batch1_s1000_inq4_forprofit.mat
% ggbs_t = zeros(20,1);
% for i = 1:length(prob_set)
%     p = prob_set{i};
%     for j = 1:size(p,2)
%         [~,id]=sort(p(:,j),'descend');
% %         if id(1)==1705&&id(2)==1704
%         if id(1) == 1705
%             ggbs_t(i) =  j;
%             break;
%         end
%     end
% end


% z = zeros(20,1);
% for i = 1:20
%     ppp = prob_set{i};
%     for j = 1:100
%         if find(ppp(:,j)==max(ppp(:,j))) == 1705
%             ppp(1705,j) = 0;
%             if find(ppp(:,j)==max(ppp(:,j))) == 1704
%                 z(i) = j;
%                 break;
%             end
%         end
%     end
% end

%% best probability
% div = zeros(20,100);
% for i = 1:20
%     ppp = prob_set{i};
%     for j = 1:100
%         div(i,j) = JSdivergence(target,ppp(:,j));
%     end
% end
meanBP = [7.37181630639441	7.29274617024863	7.18984690314466	7.05669778363843	6.94442000650883	6.80974743624641	6.66382194208975	6.51021683398286	6.35588980656919	6.20436099165340	6.01700127969016	5.85123422831940	5.68839605538334	5.40082549903438	5.17224013573265	4.99171058906081	4.80928906476675	4.68213497522559	4.51952433863063	4.38476262190629	4.29845403720706	4.23343641106740	4.16853353920547	4.13598225275117	4.10047415931933	4.07347400552143	4.06208926343539	4.04472883475088	4.02709516020765	4.02513768601080	4.00273455898411	3.98609944015119	3.97032665300535	3.97787944523077	3.95023419980989	3.93758137933846	3.91804939881442	3.92107181458831	3.89964481919931	3.89441121392157	3.89376174555952	3.86430380950268	3.87813995912419	3.86005056858307	3.86270660089557	3.83863697927057	3.83121426093956	3.83529811424524	3.82471513886434	3.82348092053100	3.81029644625408	3.80893751083120	3.79286808068381	3.79564238700985	3.79393025891387	3.78182544032359	3.76985630695501	3.77633151264410	3.77533493562407	3.76158554919837	3.76043619104898	3.75427018505585	3.74135108747457	3.74654591892668	3.74859783650145	3.74247235348386	3.72594669564970	3.72970178578540	3.73345090164738	3.72438507588340	3.73887096945811	3.72002257693624	3.71065559921310	3.71457534489914	3.69833940588053	3.69930863257704	3.68625386480637	3.70181631967748	3.69347574813450	3.67476356554126	3.68114120769721	3.68922028875841	3.66435809247925	3.66571916267219	3.66770290773158	3.66789559778746	3.66194479857713	3.66990296806033	3.65561895712120	3.64771244865367	3.65398966828031	3.65175609654388	3.64819821730311	3.64637956604549	3.64860928360260	3.63715263146683	3.63222482883375	3.64229768449333	3.63955377684072	3.63882189363281];
stdBP = [0.00938578532571113	0.0207536518029259	0.0381296771232585	0.0654828196637642	0.0721310617942162	0.100997996206036	0.168286038284641	0.223814859208657	0.238660194253662	0.271657902223360	0.301718101336298	0.336502831939488	0.353114003066352	0.411251026179252	0.405684299439039	0.387053404019995	0.409625948137384	0.431321806516870	0.425156371111287	0.444793485894508	0.438046381978787	0.438416143318199	0.435787699921875	0.398854278798832	0.380881798839248	0.377688566988973	0.377541301104963	0.386575306645437	0.370653986540270	0.372824701023028	0.394534610440128	0.389605036664722	0.397662682404012	0.379797075163332	0.387120951650180	0.380834333303616	0.372276232247898	0.381027117481754	0.396932344652836	0.378529964031373	0.381683462712829	0.375639370433908	0.372428665812460	0.381869133555625	0.382202587643885	0.367954236663633	0.387060141141238	0.377636986437527	0.378428756609270	0.370444130437420	0.369194476586549	0.367923053211939	0.381176463507250	0.370839225674228	0.373641680984173	0.371788430503543	0.390665304641283	0.359345761697469	0.377198948857140	0.372294575004108	0.374154917563112	0.366782725003600	0.371242434709967	0.370542402196341	0.381957460029583	0.359567245153770	0.377387703666508	0.365768725720236	0.363850088694115	0.373242746528084	0.363592278489050	0.355485021737495	0.357447678535265	0.364222725718696	0.375675892207855	0.373633224336192	0.364835020394739	0.380543520013846	0.366956861193018	0.379935484267150	0.374834215143025	0.372741409144341	0.376853864155613	0.375606257821061	0.386474275311576	0.372908073287333	0.367370107288296	0.364213494779684	0.384042545604097	0.375925632983470	0.373157564115282	0.369651474167298	0.370235397712232	0.376016746620416	0.382767607930213	0.375986661635555	0.362769590318639	0.365364289991145	0.372487894307964	0.375418058412212];

target = sparse(zeros(2455,1));
target(1704) = 0.3709;
target(1705) = 0.6291;

%% ggbs individual test results
%rng(6)
load ggbs_test6_s10000_inq10_forprofit.mat;
ggbs6 = zeros(1,49);
prob = prob_set{1};
for i = 1:49
    ggbs6(i) = JSdivergence(target,prob(:,i));
end
for j = 1:49
    [~,id]=sort(prob(:,j),'descend');
    if id(1)==1705
        ggbs6_t =  j;
        break;
    end
end

%rng(1)
load ggbs_test1_s10000_inq10_forprofit.mat;
ggbs1 = zeros(1,41);
prob = prob_set{1};
for i = 1:41
    ggbs1(i) = JSdivergence(target,prob(:,i));
end
for j = 1:41
    [~,id]=sort(prob(:,j),'descend');
    if id(1)==1705
        ggbs1_t =  j;
        break;
    end
end

% rng(11)
load ggbs_test11_s10000_inq10_forprofit.mat;
ggbs11 = zeros(1,54);
prob = prob_set{1};
for i = 1:54
    ggbs11(i) = JSdivergence(target,prob(:,i));
end
for j = 1:54
    [~,id]=sort(prob(:,j),'descend');
    if id(1)==1705
        ggbs11_t =  j;
        break;
    end
end

% rng(16)
load ggbs_test16_s10000_inq10_forprofit.mat;
ggbs16 = zeros(1,40);
for i = 1:40
    ggbs16(i) = JSdivergence(target,prob(:,i));
end
for j = 1:40
    [~,id]=sort(prob(:,j),'descend');
    if id(1)==1705
        ggbs16_t =  j;
        break;
    end
end

% rng(36)
load ggbs_test36_s10000_inq10_forprofit.mat;
ggbs36 = zeros(1,40);
for i = 1:40
    ggbs36(i) = JSdivergence(target,prob(:,i));
end
for j = 1:40
    [~,id]=sort(prob(:,j),'descend');
    if id(1)==1705
        ggbs36_t =  j;
        break;
    end
end

% rng(31-36)
load ggbs_batch_31_36_s1000_inq10_forprofit.mat;
ggbs3136 = zeros(6,40);
ggbs3136_t = zeros(6,1);
s = 1e4;
c = Xf(:,26:30)*price'-cv;
for j = 1:6
    prob = prob_set{j};
    pairs = pairs_set{j};
    parfor i = 1:40
        alg = 'profit';
        [A,W] = appObjDistribution(s,d,w,Xf,dX,dXID,c,[],alg,pairs(1:i,:),[]);
        probability_obj = A/sum(A);
        ggbs3136(j,i) = JSdivergence(target,probability_obj);
%         [~,id]=sort(prob(:,i),'descend');
%         if id(1)==1705
%             ggbs3136_t(j) =  i;
%             break;
%         end
    end
end

% rng(41-46)
load ggbs_batch_41_46_s1000_inq10_forprofit.mat;
ggbs4146 = zeros(6,40);
ggbs4146_t = zeros(6,1);
for j = 1:6
    prob = prob_set{j};
    for i = 1:40
        ggbs4146(j,i) = JSdivergence(target,prob(:,i));
%         [~,id]=sort(prob(:,i),'descend');
%         if id(1)==1705
%             ggbs4146_t(j) =  i;
%             break;
%         end
    end
end

% rng(51-53)
load ggbs_batch_51_53_s1000_inq10_forprofit.mat;
ggbs5153 = zeros(3,NQ);
ggbs5153_t = zeros(3,1);
for j = 1:3
    prob = prob_set{j};
    for i = 2:(NQ+1)
        ggbs5153(j,i-1) = JSdivergence(target,prob(:,i));
%         [~,id]=sort(prob(:,i),'descend');
%         if id(1)==1705
%             ggbs5153_t(j) =  i;
%             break;
%         end
    end
end

%% pref individual test results
% rng(0)
load pref_test0_s1000_inq10_forpref_only.mat;
pref0 = zeros(1,40);
s = 1e4;
c = Xf(:,26:30)*price'-cv;
parfor i = 1:40
    alg = 'profit';
    [A,W] = appObjDistribution(s,d,w,Xf,dX,dXID,c,[],alg,pairs(1:i,:),[]);
    probability_obj = A/sum(A);
    pref0(i) = JSdivergence(target,probability_obj);
end
for j = 1:40
    [~,id]=sort(prob(:,j),'descend');
    if id(1)==1705
        pref0_t =  j;
        break;
    end
end

% rng(4)
load pref_test4_s1000_inq10_forpref_only.mat;
pref4 = zeros(1,40);
parfor i = 1:40
    alg = 'profit';
    [A,W] = appObjDistribution(s,d,w,Xf,dX,dXID,c,[],alg,pairs(1:i,:),[]);
    probability_obj = A/sum(A);
    pref4(i) = JSdivergence(target,probability_obj);
end
for j = 1:40
    [~,id]=sort(prob(:,j),'descend');
    if id(1)==1705
        pref4_t =  j;
        break;
    end
end

% rng(5)
load pref_test5_s1000_inq10_forpref_only.mat;
pref5 = zeros(1,40);
parfor i = 1:40
    alg = 'profit';
    [A,W] = appObjDistribution(s,d,w,Xf,dX,dXID,c,[],alg,pairs(1:i,:),[]);
    probability_obj = A/sum(A);
    pref5(i) = JSdivergence(target,probability_obj);
end
for j = 1:40
    [~,id]=sort(prob(:,j),'descend');
    if id(1)==1705
        pref5_t =  j;
        break;
    end
end
% rng(51-68)
load pref_test_51_68_s1000_inq10_forpref_only.mat;
pref5168 = zeros(17,NQ);
pref5168_t = zeros(17,1);
for j = 1:17
    prob = prob_set{j};
    for i = 2:(NQ+1)
        pref5168(j,i-1) = JSdivergence(target,prob(:,i));
%         [~,id]=sort(prob(:,i),'descend');
%         if id(1)==1705
%             pref5168_t(j) =  i;
%             break;
%         end
    end
end

% from Toubia for pref only
load pref_simple_batch_1_20_s10000_inq10_forpref_only.mat
pref_simple = zeros(1,40);
ppp = prob_set{1};
parfor i = 1:40
    alg = 'profit';
    pref_simple(i) = JSdivergence(target,ppp(:,i));
end

ggbs_div = [ggbs1(2:31);ggbs6(2:31);ggbs11(2:31);...
    ggbs16(2:31);ggbs36(2:31)];
%     ggbs3136(:,2:31);ggbs4146(:,2:31);...
%     ggbs5153];

pref_div = [pref0;pref4;pref5];
pref_div = [pref_div(:,1:30);pref5168];
figure;hold on;
shadedErrorBar(1:30,meanBP(2:31),stdBP(2:31),'k',1); 
shadedErrorBar(1:30,pref_div,{@mean,@std},'-sk',1); 
shadedErrorBar(1:30,ggbs_div,{@mean,@std},'k',1); 
plot(meanBP(2:31),'-ok','LineWidth',3);
% plot(mean(pref_div),'-sk','LineWidth',3);
plot(pref_simple(2:31),'--','LineWidth',3);
plot(mean(ggbs_div),'-k','LineWidth',3);


%% check w convergence
% calculate w under all information
A = bsxfun(@times,dX,2*(dX*w'>0)-1);
A(abs(dX*w')<1e-3,:) = [];
addpath('C:\Users\Max Yi Ren\Documents\Research\Code\Tools\liblinear\matlab');
model = train([ones(size(A,1),1);-ones(size(A,1),1)],...
                    sparse([A;-A]),'-s 0 -c 1e6 -e 1e-6 -q');
w_all = model.w;

% calculate w under best probability approach
load bestprob_batch0_s10000_inq100_forprofit.mat
dist_w_bestprob = zeros(20,99);
DX = dXID + dXID';
for i = 1:20
    pairi = pairs_set{i};
    for j = 1:99
        Ai = dX(DX(sub2ind(size(DX),pairi(1:j,1),pairi(1:j,2))),:);
        Ai = bsxfun(@times,Ai,sign(pairi(1:j,2)-pairi(1:j,1)));
        modeli = train([ones(size(Ai,1),1);-ones(size(Ai,1),1)],...
                        sparse([Ai;-Ai]),'-s 0 -c 1e6 -e 1e-6 -q');
        wi = modeli.w;
        dist_w_bestprob(i,j) = wi/norm(wi)*w'/norm(w);
    end
end

% calculate w under ggbs approach
load ggbs_test1_s10000_inq10_forprofit.mat
dist_w_ggbs1 = zeros(1,NQ);
pairi = pairs_set{1};
for j = 1:NQ
    Ai = dX(DX(sub2ind(size(DX),pairi(1:j,1),pairi(1:j,2))),:);
    Ai = bsxfun(@times,Ai,sign(pairi(1:j,2)-pairi(1:j,1)));
    modeli = train([ones(size(Ai,1),1);-ones(size(Ai,1),1)],...
                    sparse([Ai;-Ai]),'-s 0 -c 1e6 -e 1e-6 -q');
    wi = modeli.w;
    dist_w_ggbs1(j) = wi/norm(wi)*w'/norm(w);
end
load ggbs_test6_s10000_inq10_forprofit.mat
dist_w_ggbs6 = zeros(1,NQ);
pairi = pairs_set{1};
for j = 1:NQ
    Ai = dX(DX(sub2ind(size(DX),pairi(1:j,1),pairi(1:j,2))),:);
    Ai = bsxfun(@times,Ai,sign(pairi(1:j,2)-pairi(1:j,1)));
    modeli = train([ones(size(Ai,1),1);-ones(size(Ai,1),1)],...
                    sparse([Ai;-Ai]),'-s 0 -c 1e6 -e 1e-6 -q');
    wi = modeli.w;
    dist_w_ggbs6(j) = wi/norm(wi)*w'/norm(w);
end
load ggbs_test11_s10000_inq10_forprofit.mat
dist_w_ggbs11 = zeros(1,NQ);
pairi = pairs_set{1};
for j = 1:NQ
    Ai = dX(DX(sub2ind(size(DX),pairi(1:j,1),pairi(1:j,2))),:);
    Ai = bsxfun(@times,Ai,sign(pairi(1:j,2)-pairi(1:j,1)));
    modeli = train([ones(size(Ai,1),1);-ones(size(Ai,1),1)],...
                    sparse([Ai;-Ai]),'-s 0 -c 1e6 -e 1e-6 -q');
    wi = modeli.w;
    dist_w_ggbs11(j) = wi/norm(wi)*w'/norm(w);
end
load ggbs_test16_s10000_inq10_forprofit.mat
dist_w_ggbs16 = zeros(1,NQ);
pairi = pairs;
for j = 1:NQ
    Ai = dX(DX(sub2ind(size(DX),pairi(1:j,1),pairi(1:j,2))),:);
    Ai = bsxfun(@times,Ai,sign(pairi(1:j,2)-pairi(1:j,1)));
    modeli = train([ones(size(Ai,1),1);-ones(size(Ai,1),1)],...
                    sparse([Ai;-Ai]),'-s 0 -c 1e6 -e 1e-6 -q');
    wi = modeli.w;
    dist_w_ggbs16(j) = wi/norm(wi)*w'/norm(w);
end
load ggbs_test36_s10000_inq10_forprofit.mat
dist_w_ggbs36 = zeros(1,NQ);
pairi = pairs;
for j = 1:NQ
    Ai = dX(DX(sub2ind(size(DX),pairi(1:j,1),pairi(1:j,2))),:);
    Ai = bsxfun(@times,Ai,sign(pairi(1:j,2)-pairi(1:j,1)));
    modeli = train([ones(size(Ai,1),1);-ones(size(Ai,1),1)],...
                    sparse([Ai;-Ai]),'-s 0 -c 1e6 -e 1e-6 -q');
    wi = modeli.w;
    dist_w_ggbs36(j) = wi/norm(wi)*w'/norm(w);
end
load ggbs_batch_31_36_s1000_inq10_forprofit.mat
dist_w_ggbs3136 = zeros(6,NQ);
for i = 1:6
    for j = 1:NQ
        Ai = dX(DX(sub2ind(size(DX),pairi(1:j,1),pairi(1:j,2))),:);
        Ai = bsxfun(@times,Ai,sign(pairi(1:j,2)-pairi(1:j,1)));
        modeli = train([ones(size(Ai,1),1);-ones(size(Ai,1),1)],...
                        sparse([Ai;-Ai]),'-s 0 -c 1e6 -e 1e-6 -q');
        wi = modeli.w;
        dist_w_ggbs3136(i,j) = wi/norm(wi)*w'/norm(w);
    end
end
load ggbs_batch_41_46_s1000_inq10_forprofit.mat
dist_w_ggbs4146 = zeros(6,NQ);
for i = 1:6
    for j = 1:NQ
        Ai = dX(DX(sub2ind(size(DX),pairi(1:j,1),pairi(1:j,2))),:);
        Ai = bsxfun(@times,Ai,sign(pairi(1:j,2)-pairi(1:j,1)));
        modeli = train([ones(size(Ai,1),1);-ones(size(Ai,1),1)],...
                        sparse([Ai;-Ai]),'-s 0 -c 1e6 -e 1e-6 -q');
        wi = modeli.w;
        dist_w_ggbs4146(i,j) = wi/norm(wi)*w'/norm(w);
    end
end
dist_w_ggbs = [dist_w_ggbs1;dist_w_ggbs6;dist_w_ggbs11;...
    dist_w_ggbs16;dist_w_ggbs36;];
%     dist_w_ggbs3136;dist_w_ggbs4146];

% calculate w under prefonly approach
load pref_test0_s1000_inq10_forpref_only.mat
dist_w_pref0 = zeros(1,NQ);
pairi = pairs;
for j = 1:NQ
    Ai = dX(DX(sub2ind(size(DX),pairi(1:j,1),pairi(1:j,2))),:);
    Ai = bsxfun(@times,Ai,sign(pairi(1:j,2)-pairi(1:j,1)));
    modeli = train([ones(size(Ai,1),1);-ones(size(Ai,1),1)],...
                    sparse([Ai;-Ai]),'-s 0 -c 1e6 -e 1e-6 -q');
    wi = modeli.w;
    dist_w_pref0(j) = wi/norm(wi)*w'/norm(w);
end
load pref_test4_s1000_inq10_forpref_only.mat
dist_w_pref4 = zeros(1,NQ);
pairi = pairs;
for j = 1:NQ
    Ai = dX(DX(sub2ind(size(DX),pairi(1:j,1),pairi(1:j,2))),:);
    Ai = bsxfun(@times,Ai,sign(pairi(1:j,2)-pairi(1:j,1)));
    modeli = train([ones(size(Ai,1),1);-ones(size(Ai,1),1)],...
                    sparse([Ai;-Ai]),'-s 0 -c 1e6 -e 1e-6 -q');
    wi = modeli.w;
    dist_w_pref4(j) = wi/norm(wi)*w'/norm(w);
end
load pref_test5_s1000_inq10_forpref_only.mat
dist_w_pref5 = zeros(1,NQ);
pairi = pairs;
for j = 1:NQ
    Ai = dX(DX(sub2ind(size(DX),pairi(1:j,1),pairi(1:j,2))),:);
    Ai = bsxfun(@times,Ai,sign(pairi(1:j,2)-pairi(1:j,1)));
    modeli = train([ones(size(Ai,1),1);-ones(size(Ai,1),1)],...
                    sparse([Ai;-Ai]),'-s 0 -c 1e6 -e 1e-6 -q');
    wi = modeli.w;
    dist_w_pref5(j) = wi/norm(wi)*w'/norm(w);
end
load pref_test_51_68_s1000_inq10_forpref_only.mat
dist_w_pref5168 = zeros(17,NQ);
for i = 1:17
    pairi = pairs_set{i};
    for j = 1:NQ
        Ai = dX(DX(sub2ind(size(DX),pairi(1:j,1),pairi(1:j,2))),:);
        Ai = bsxfun(@times,Ai,sign(pairi(1:j,2)-pairi(1:j,1)));
        modeli = train([ones(size(Ai,1),1);-ones(size(Ai,1),1)],...
                        sparse([Ai;-Ai]),'-s 0 -c 1e6 -e 1e-6 -q');
        wi = modeli.w;
        dist_w_pref5168(i,j) = wi/norm(wi)*w'/norm(w);
    end
end

load pref_simple_batch_1_20_s10000_inq10_forpref_only.mat
dist_w_pref_simple = zeros(1,NQ);
pairi = pairs_set{1};
for j = 1:NQ
    Ai = dX(DX(sub2ind(size(DX),pairi(1:j,1),pairi(1:j,2))),:);
    Ai = bsxfun(@times,Ai,sign(pairi(1:j,2)-pairi(1:j,1)));
    modeli = train([ones(size(Ai,1),1);-ones(size(Ai,1),1)],...
                    sparse([Ai;-Ai]),'-s 0 -c 1e6 -e 1e-6 -q');
    wi = modeli.w;
    dist_w_pref_simple(j) = wi/norm(wi)*w'/norm(w);
end

figure;hold on;
shadedErrorBar(1:30,dist_w_bestprob(:,1:NQ),{@mean,@std},'k',1); 
shadedErrorBar(1:30,[dist_w_pref0;dist_w_pref4;dist_w_pref5],{@mean,@std},'k',1);
shadedErrorBar(1:30,dist_w_ggbs,{@mean,@std},'k',1);

% plot(mean(dist_w_bestprob(:,1:NQ)));
% plot(mean([dist_w_pref0;dist_w_pref4;dist_w_pref5]));
plot(dist_w_pref_simple);
plot(mean(dist_w_ggbs));