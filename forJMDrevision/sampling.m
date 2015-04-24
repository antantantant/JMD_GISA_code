function [W,w0,C,weights] = sampling(s,d,W0,A)
    
    if isempty(A)
        W = mvnrnd(zeros(1,d),eye(d),s); % generate candidates
        w0 = zeros(1,d);
        C = 1;
        weights = 1;
    else
        [w0, C, weights] = crossvalidation(sparse(A),ones(size(A,1),1)); 
        
        step = length(C);
        W = zeros(size(W0,1)*step,size(W0,2));
        for i = 1:step
            As = bsxfun(@times, A, sqrt(exp(A*w0(i,:)'))./(1+exp(A*w0(i,:)')));
            As(isnan(As))=0;
            Sigma_inv = (eye(d)/C(i)+(As'*As));
%             Sigma_inv = (eye(d)/C(i)+(A'*A));
            cond(A'*A)
            W_ = W0*(inv(chol(Sigma_inv)))';
            W_ = bsxfun(@plus, W_, w0(i,:));  
            W((i-1)*size(W0,1)+(1:size(W0,1)),:) = W_;        
        end
        
%         step = 10;
%         W = zeros(size(W0,1)*step,size(W0,2));
%         W_base = zeros(size(W0,1),size(W0,2),length(C));
%         for i = 1:length(C)
%             As = bsxfun(@times, A, sqrt(exp(A*w0(i,:)'))./(1+exp(A*w0(i,:)')));
%             As(isnan(As))=0;
%             Sigma_inv = (eye(d)/C(i)+(As'*As));
%             W_ = W0*(inv(chol(Sigma_inv)))';
%             W_ = bsxfun(@plus, W_, w0(i,:));  
%             W_base(:,:,i) = W_;        
%         end
%         
%         if length(C)==2
%             for i = 1:step
%                 a = (step-i)/(step-1);
%                 W((i-1)*size(W0,1)+(1:size(W0,1)),:) = a*W_base(:,:,1)+(1-a)*W_base(:,:,2);         
%             end
%         else
%             W = W_base(:,:,1);
%         end

    end