function xoverKids  = gacrossover(parents,options,GenomeLength,FitnessFcn,unused,thisPopulation)
%CROSSOVERSCATTERED Position independent crossover function.
%   XOVERKIDS = CROSSOVERSCATTERED(PARENTS,OPTIONS,GENOMELENGTH, ...
%   FITNESSFCN,SCORES,THISPOPULATION) creates the crossover children XOVERKIDS
%   of the given population THISPOPULATION using the available PARENTS.
%   In single or double point crossover, genomes that are near each other tend
%   to survive together, whereas genomes that are far apart tend to be
%   separated. The technique used here eliminates that effect. Each gene has an
%   equal chance of coming from either parent. This is sometimes called uniform
%   or random crossover.
%
%   Example:
%    Create an options structure using CROSSOVERSCATTERED as the crossover
%    function
%     options = gaoptimset('CrossoverFcn' ,@crossoverscattered);

%   Copyright 2003-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2012/08/21 00:21:36 $


% How many children to produce?
nKids = length(parents)/2;
% Allocate space for the kids
xoverKids = zeros(nKids,GenomeLength);

% To move through the parents twice as fast as thekids are
% being produced, a separate index for the parents is needed
index = 1;
% for each kid...
for i=1:nKids
    % get parents
    r1 = parents(index);
    index = index + 1;
    r2 = parents(index);
    index = index + 1;
    % Randomly select half of the genes from each parent
    % This loop may seem like brute force, but it is twice as fast as the
    % vectorized version, because it does no allocation.
    for j = 1:6
        if(rand > 0.5)
            xoverKids(i,(j-1)*5+(1:5)) = thisPopulation(r1,(j-1)*5+(1:5));
        else
            xoverKids(i,(j-1)*5+(1:5)) = thisPopulation(r2,(j-1)*5+(1:5));
        end
    end
end
