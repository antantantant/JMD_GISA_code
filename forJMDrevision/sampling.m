function [W,w0] = sampling(s,d,W0,Xf,dX,dXID,pairs,C)
    
    if isempty(pairs)
        W = mvnrnd(zeros(1,d),eye(d),s); % generate candidates
        w0 = zeros(1,d);
    else
        theta = 1; % logistic penalty
        A = dX(dXID(sub2ind(size(dXID),pairs(:,1),pairs(:,2))),:);
        A = bsxfun(@times,A,sign(pairs(:,2)-pairs(:,1)));
        b = zeros(size(pairs,1),1);
            
        model = train(ones(size(A,1),1),...
                            sparse(theta*A), sprintf('-s 0 -e 1e-6 -q -c %f', C));
        w0 = model.w;
%         w00 = ((A'*A+1/C*eye(size(A,2)))\A'*ones(size(A,1),1))';
%             [w0,C] = crossvalidation(sparse(A),ones(size(A,1),1));
%         W = myMetropolisHasting(w0,s,s,A,C,theta);

        As = bsxfun(@times, theta*A, sqrt(exp(theta*A*w0'))./(1+exp(theta*A*w0')));
        Sigma_inv = (eye(d)/C+(As'*As));
        W = W0*(inv(chol(Sigma_inv)))';
        W = bsxfun(@plus, W, w0);
        
%             W = mvnrnd(w0,inv(Sigma),s);
% %             W = feasible(W,A,b);
% %             if isempty(W)
% %                 stop = 1;
% %             end

    end
        
    function WC = feasible(WC,A,b)
        feasible = sum(bsxfun(@lt, A*WC', b),1)==0;
        WC = WC(feasible,:);