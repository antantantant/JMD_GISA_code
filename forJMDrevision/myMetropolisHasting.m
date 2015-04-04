%my metropolis hasting algorithm

function [smpl,accept] = myMetropolisHasting(start,nsamples,burnin,A,C,theta)

    % Put the replicates dimension second.
    distnDims = size(start,2);
    smpl = zeros([nsamples,distnDims],'double');

    x0 = start;  %x0  is the place holder for the current value
    accept =0;    
    % Metropolis-Hasting Algorithm.
    U = log(rand(1,nsamples+burnin));
    for i = 1-burnin:nsamples
        y = proprnd(x0); % sample from proposal dist'n
        %save the evaluation time for symmetric proposal dist'n
        rho = logpdf(y',A,C,theta)-logpdf(x0',A,C,theta); 
        % Accept or reject the proposal.
        Ui = U(:,i+burnin);
        acc = (Ui<= min(rho,0)) + 0;
        if acc>0
            x0(acc,:) = y(acc,:); % preserves x's shape.
        end
        accept = accept + acc(1);
        if i>0 % burnin
            smpl(i,:) = x0;
        end;
    end;

    % Accept rate can be used to optimize the choice of scale parameters in
    % random walk MH sampler. See for example Roberts, Gelman and Gilks (1997).
    accept = accept/(nsamples+burnin);