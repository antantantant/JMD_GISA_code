function x_new = proprnd(x, An, dn, C)
global eee
p = size(An,2);
n = size(An,1);
eee = chi2rnd(ones(n,1));
model = train(eee, sparse(ones(n,1)), sparse(An), ...
        sprintf('-s 0 -e 1e-6 -q -c %f', C));
x_new = (model.w);
% x_new = ((An'*An+1e-3*eye(p))\An'*e)';

%     x_new = mvnrnd(w0, Sigma, 1);  