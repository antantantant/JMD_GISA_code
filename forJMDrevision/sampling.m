function [W,w0,C,weights] = sampling(s,d,A)
    
    if isempty(A)
        W = mvnrnd(zeros(1,d),eye(d),s); % generate candidates
        w0 = zeros(1,d);
        C = 1;
        weights = ones(s,1);
    else
        % update estimation of w, row vector
        [w0, C] = crossvalidation(sparse(A),ones(size(A,1),1)); 
        % calculate sigma
        As = bsxfun(@times, A, sqrt(exp(A*w0'))./(1+exp(A*w0')));
        As(isnan(As))=0;
        Hessian = (eye(length(w0))/C+(As'*As));
        
        
         %% importance sampling
%         % first figure out the constant
% %         const = 1;
% %         Sigmac = Sigma;
% %         s0 = 1e3;
% %         pqW = 1;
%         
%         
% %         f = @(x)bisectfunc(x, w0, Sigma, A, C);
%         a = 1; b = 1e-32;
%         if (bisectfunc(a,w0,Sigma,A, C) * bisectfunc(b,w0,Sigma,A, C))>0 
%             disp('Wrong choice bro')
%         else
%             p = (a + b)/2;
%             err = abs(bisectfunc(p, w0, Sigma, A, C));
%             while err > 1e-3
%                 if bisectfunc(a, w0, Sigma, A, C)*bisectfunc(p, w0, Sigma, A, C)<0 
%                    b = p;
%                 else
%                    a = p;          
%                 end
%                 p = (a + b)/2; 
%                 err = abs(bisectfunc(p, w0, Sigma, A, C));
%             end
%         end
%         const = p;
% %         while abs(pqW)>1e-3
% % %             W = mvnrnd(w0,Sigmac,s0);
% % %             qW = log(mvnpdf(W, w0, Sigmac));
% % %             dW = sum(W.^2,2);
% % %             temp = -W*A';
% % %             temp2 = log(1+exp(temp));
% % %             temp2(temp2==inf) = temp(temp2==inf);
% % %             pW = -dW/2/C+sum(-temp2,2);
% % %             pqW = pW-qW;
% % %             y = max(pqW);
% % %             const = y+log(sum(exp(bsxfun(@minus,pqW,y))))-log(s0);
% % 
% %             if pqW<-1e-3
% %                 break;
% %             elseif pqW>1e-3
% %                 const = const/2;
% %             end
% % %             const = pqW;
% % %             Sigmac = Sigmac*exp(const);
% %         end
%         
%         Sigmac = Sigma*const;
%         W = mvnrnd(w0,Sigmac,s);
%         qW = log(mvnpdf(W, w0, Sigmac));
%         dW = sum(W.^2,2);
%         temp = -W*A';
%         temp2 = log(1+exp(temp));
%         temp2(temp2==inf) = temp(temp2==inf);
%         pW = -dW/2/C+sum(-temp2,2)-log(const);
%         pqW = pW-qW;
%         weights = exp(pqW);
%         if sum(weights)==0
%             weights=ones(s,1);
%         end
%         weights = weights/sum(weights);

        %% MCMC
%         c = 1e1;
%         Sigma = inv(Hessian);
%         Sigma = Sigma*1e-3+w0'*w0*c; 
% 
%         [V,D] = eig(Sigma);
%         ll = max(diag(D));
%         ss = sqrt(ll);
%         dd = 3*w0/norm(w0)*ss;
%         logproppdf(mvnrnd(w0, Sigma, 1),[],w0,Sigma)-logproppdf(w0,[],w0,Sigma)
%         logpdf(mvnrnd(w0, Sigma, 1),A,C)-logpdf(w0,A,C)
        
%         accept = 0;
%         s_t = 1e3;
%         c = 0.01;
%         co = c*Sigma;
%         count = 0; limit = 10;
%         while count < limit && (accept>0.4 || accept<0.2)
%             [W, accept] = mhsample(w0, s_t, 'logpdf',@(x)logpdf(x,A,C), ...
%                 'proprnd',@(x)proprnd(x,w0,co), ...
%                 'logproppdf',@(x,y)logproppdf(x,y,w0,co));
%             if accept>0.4
%                 co = co*100;
%             elseif accept<0.2
%                 co = co/10;
%             end
%             count = count + 1;
%         end
%         accept


        p = numel(w0);
        if size(A,1)>1
            AA = (licols(A'))';
        else
            AA = A;
        end
        d = AA*w0';
        [~,st] = sort(d);
        d_st = d(st);
        AAA = AA(st,:);
%         An = AAA(1:min(size(AAA,1),p),:);
%         dn = d_st(1:min(size(AAA,1),p));
        An = AAA;
        dn = d_st;
        [W, accept] = mhsample(w0, s, 'logpdf',@(x)logpdf(x,A,C), ...
                'proprnd',@(x)proprnd(x,An,dn,C), ...
                'logproppdf',@(x,y)logproppdf(x,y,An,dn));
        accept
% 
%         [W, accept] = mhsample(w0, s, 'logpdf',@(x)logpdf(x,A,C), ...
%                 'proprnd',@(x)proprnd(x,w0,Sigma), ...
%                 'logproppdf',@(x,y)logproppdf(x,y,w0,Sigma));
%         accept

%         Sigma = inv(Hessian);
%         Sigma = Sigma+w0'*w0*1e4; 
%         Sigma = Sigma/max(diag(Sigma));
%         OPT.logpdf = @(x)logpdf(x',A,C);    % Handle to log-density function
%         OPT.V = Sigma;             % Covariance matrix
%         OPT.Dims = d;        % Number of dimensions
%         OPT.dsp = 1;            % Display status of MCMC runs, default value = 1
%         OPT.Mmax=1e4;           % Number of MCMC samples, recommend 1e6 sample for slice sampling
%         OPT.div = 1;
%         OPT.hgrd= 1e2;          % Number of grid point in histograms, must be multiple of 10
%         OPT.filename='AM_output';
%         W = AM(w0',OPT)';
%         W1 = AM(rand(d,1),OPT);
        weights=ones(s,1);
       
    end