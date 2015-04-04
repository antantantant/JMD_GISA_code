% test metropolis hasting
close;
figure;hold on;
A = [1,1];
C = 1;
W = myMetropolisHasting_mex([0,0],1e4,1e4,A,C);
x = 0:0.1:10;
y = 0:0.1:10;
[X,Y] = meshgrid(x,y);
Z = zeros(length(x),length(y));
for i = 1:length(x)
    for j = 1:length(y)
        Z(i,j) = exp(logpdf([x(i);y(j)],A,C));
    end
end
contourf(X,Y,Z);
scatter(W(:,1),W(:,2));