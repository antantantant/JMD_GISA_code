
% mex mvtrcnml.c
rng(0);
s = 1e3;
p = 2;
A = [1,-10;-1,20;1,-1;1,-2];
b = [0;0;0;0];
w0 = [1,0.7];
tic
W = mvtrcnml(s,p,A',b,w0);
W = W(:,s+1:end)';
toc
figure;
scatter(W(:,1),W(:,2));

tic
W2 = mvtnconstrained(s,p,[],[],A,b,w0);
toc