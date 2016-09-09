% creat artificial response data in 2d
n_repeat = 100;
A = [1,-0.6;-0.6,1];
A = repmat(A,n_repeat,1);
% training
[w0, C] = crossvalidation(sparse(A),ones(size(A,1),1)); 

%visualize logprob and w0
maxx = 500*round(w0(1));
maxy = 500*round(w0(2));
intx = (maxx+1)/100;
inty = (maxy+1)/100;
[X,Y] = meshgrid(-1:intx:maxx, -1:inty:maxy);
Z = zeros(size(X,1),size(Y,2));
for i = 1:size(X,1)
    for j = 1:size(Y,2)
        w = [X(i,j),Y(i,j)];
        Z(i,j) = exp(logpdf(w,A,C));
    end
end
figure;hold on;
surf(X,Y,Z);xlabel('x');ylabel('y');
% contourf(X,Y,Z);xlabel('x');ylabel('y');
plot(w0(1),w0(2),'.','MarkerSize',40);

% try mcmc idea
s = 1e4;
p = numel(w0);
A = (licols(A'))';
d = A*w0';
[~,st] = sort(d);
d_st = d(st);
A = A(st,:);
An = A(1:max(size(A,1),p),:);
dn = d_st(1:max(size(A,1),p));
[W, accept] = mhsample(w0, s, 'logpdf',@(x)logpdf(x,A,C), ...
        'proprnd',@(x)proprnd(x,An,dn), ...
        'logproppdf',@(x,y)logproppdf(x,y,An,dn));
accept

figure; plot(W(:,1),W(:,2),'.','MarkerSize',20);

% As = bsxfun(@times, A, sqrt(exp(A*w0'))./(1+exp(A*w0')));
% As(isnan(As))=0;
% Hessian = (eye(length(w0))/C+(As'*As));
% Sigma = inv(Hessian);
