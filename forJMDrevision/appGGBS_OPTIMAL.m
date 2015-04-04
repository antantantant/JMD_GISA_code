function [probability_obj_set,pairs,method] = ...
    appGGBS_OPTIMAL(s,d,Xf,dX,dXID,c,inq,const,alg,queryID,pairs,w,target,nq)
    nx = size(Xf,1);
    probability_obj_set = [];
    nquery = sum(sum(queryID>0));
    method = zeros(10001,1);
    gaoptions = gaoptimset('CreationFcn',@gacreation,...
        'CrossoverFcn',@gacrossover,'MutationFcn',@gamutate,...
        'PopulationSize',inq,'PopulationType','custom',...
        'PopInitRange',[-1;1],'Generations',100,'StallGenLimit',50);
    % find the next query
    while  nq>=0
        % calculate A
        [A,W,I,unique_I] = appObjDistribution(s,d,w,Xf,dX,dXID,c,const,alg,pairs,[]);
        probability_obj = A/sum(A);
        fprintf('%d, %.2f, %.2f \n',nq,JSdivergence(probability_obj,target),...
            probability_obj(1704));
        
%         if ~isempty(pairs)&&sum(probability_obj<(1/nx*1e-3))==nx-1
        if nq>1000
            probability_obj_set = [probability_obj_set, probability_obj];
            fprintf('max query number exceeded. terminated.');
            return;
        else
            current_dX = ga(@(dx)queryObj(W,I,unique_I,dx',nx,probability_obj),d,...
                [],[],[],[],[],[],[],gaoptions);

            u = current_dX*w;
            if u>0
                pairs = [pairs;current_dX];
            elseif u<0
                pairs = [pairs;-current_dX];
            end
            probability_obj_set = [probability_obj_set, probability_obj];
            nq = nq+1;
            nquery = nquery - 1;
        end
    end