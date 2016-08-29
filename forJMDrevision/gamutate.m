function mutationChildren = gamutate(parents,options,GenomeLength, ...
     FitnessFcn,state,thisScore,thisPopulation,scale,shrink)

     mutationChildren = zeros(length(parents),GenomeLength);
     for i=1:length(parents)
         parent = thisPopulation(parents(i),:);
         mutationChildren(i,:) = parent;
         if rand >= state.Generation/options.Generations
            ind = ceil(rand*6);
            temp = rand(1,5);
            temp_max = bsxfun(@eq,temp,max(temp,[],2));
            temp_min = bsxfun(@eq,temp,min(temp,[],2));
            temp_all = temp_max-temp_min;
            mutationChildren(i,(ind-1)*5+(1:5)) = temp_all;
         end
     end