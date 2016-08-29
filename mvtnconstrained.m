% Gibbs sampler for truncated multivariate normal distribution
% with linear constraints
% Implementation follows "Efficient Gibbs Sampling of Truncated 
% Multivariate Normal with Application to Constrained Linear Regression"
% Rodriguez-Yam, Davis, and Sharf

% Inputs
% s - sampling size
% p - dimensions
% mu - mean
% sigma - covariance
% A,b - Ax>=b
% w0 - initial feasible sampler

% current implementation for mu = zeros(p,1),
% sigma = eye(p) and b = 0
function W = mvtnconstrained(s,p,mu,sigma,A,b,w0)
    N = 2*s;
    samples = zeros(N,p);
    samples(1,:) = w0;
    for i = 2:N
%         fprintf('%d|',round(i/N*100));
        for j = 1:p
%             fprintf('.');
            AID = abs(A(:,j))>1e-3;
            if ~isempty(b(AID))
                bb = b - A(:,[1:j-1,j+1:end])*[samples(i,1:j-1),...
                    samples(i-1,j+1:end)]';
                bb = bb./A(:,j);
                bb = bb(AID);
                lb = max(bb(A(AID,j)>0));
                ub = min(bb(A(AID,j)<0));
                if isempty(lb)
                    lb = -1e2;
                end
                if isempty(ub)
                    ub = 1e2;
                end
                if abs(ub-lb)>1e-3
                    try
                        a = truncnormrnd(1,0,1,lb,ub);
%                         a = (lb+ub)/2;
                    catch
                        aa = 1;
                    end
                else
                    a = samples(i-1,j);
                end
                if isinf(a)
                    a = (lb+ub)/2;
                end
                samples(i,j) = a;
            else
                samples(i,j) = randn();
            end
        end
    end
    W = samples(s+1:end,:);
    fprintf('.');

