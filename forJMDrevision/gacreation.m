function Population = gacreation(GenomeLength,FitnessFcn,options)

    totalpopulation = sum(options.PopulationSize);
    temp = rand(totalpopulation*6,5);
    temp_max = bsxfun(@eq,temp,max(temp,[],2));
    temp_min = bsxfun(@eq,temp,min(temp,[],2));
    temp_all = temp_max-temp_min;
    Population = reshape(temp_all',30,totalpopulation)';