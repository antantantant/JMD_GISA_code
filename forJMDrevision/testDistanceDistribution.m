% mu = rand(1,24);
% R = chol(Sigma_inv);
% D = rand(24);
% R = chol(D'*D);
% sigma = eig(inv(Sigma_inv))';
rng default  % For reproducibility
% r = mvnrnd(zeros(1,24),10.^(-12:11), 1e6);
r = mvnrnd(zeros(1,30),eye(30), 1e6);

% r = r*inv(R)';
% mu = w0;
% sigma = eig(inv(Sigma_inv))'*1e5;
% rng default  % For reproducibility
% r = mvnrnd(mu,sigma,1e6);
norm(mean(r))
% D2 = bsxfun(@minus, r, mu);
D2 = sqrt(sum(r.^2,2));
figure;hist(D2);