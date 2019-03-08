%Written by Seyedali Mirjalili 

function [childGene] = mutation_continious(childGene, Pm, varMin, varMax)

[~,Gene_no] = size(childGene);

for k = 1:Gene_no
    
    R = rand();
    
    if R < Pm
        
        childGene(k) = (varMax(k) - varMin(k)) * rand() + varMin(k);
    end
end

end