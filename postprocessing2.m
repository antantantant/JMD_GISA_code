% postprocessing on when the best design is found
probs = sparse(zeros(2455,1));
load ggbs_test6_s10000_inq10_forprofit.mat;
probs = probs + prob(:,31);
load ggbs_test1_s10000_inq10_forprofit.mat;
probs = probs + prob(:,31);
load ggbs_test11_s10000_inq10_forprofit.mat;
probs = probs + prob(:,31);
load ggbs_test16_s10000_inq10_forprofit.mat;
probs = probs + prob(:,31);
load ggbs_test36_s10000_inq10_forprofit.mat;
probs = probs + prob(:,31);
load ggbs_batch_31_36_s1000_inq10_forprofit.mat;
for i = 1:6
    prob = prob_set{i};
    probs = probs + prob(:,31);
end
load ggbs_batch_41_46_s1000_inq10_forprofit.mat;
for i = 1:6
    prob = prob_set{i};
    probs = probs + prob(:,31);
end
figure;
plot(probs)


