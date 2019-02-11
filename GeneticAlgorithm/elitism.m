%Written by Seyedali Mirjalili 

function [newPopFitness,newPopGene] = elitism(nPop, popFitness, popGene, newPopFitness, newPopGene, Er)
 
Elite_no = round(nPop * Er);

[~,indx] = sort(popFitness, 'ascend');

% Worst in population replaced by elite
for k = 1 : Elite_no
    
    newPopGene(k,:) = popGene(indx(k),:);
    newPopFitness(k)  = popFitness(indx(k));
end

end