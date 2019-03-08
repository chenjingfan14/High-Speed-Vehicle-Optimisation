%Written by Seyedali Mirjalili 

function [p1Gene,p2Gene] = selection(nPop,fitness,chromosome)

% Inverting for minimisation problems
invFitness = max(fitness) - fitness;

normFitness = invFitness ./ sum(invFitness);

[~,sorted_idx] = sort(normFitness , 'descend');

tempGene = chromosome(sorted_idx,:);
tempNormFitness = normFitness(sorted_idx);

cumsum = zeros(1 , nPop);

for i = 1 : nPop
    
    for j = i : nPop
        
        cumsum(i) = cumsum(i) + tempNormFitness(j);
    end
end

R = rand(); % in [0,1]
parent1_idx = nPop;

for i = 1: nPop
    
    if R > cumsum(i)
        
        parent1_idx = i - 1;
        break;
    end
end

parent2_idx = parent1_idx;
while_loop_stop = 0; % to break the while loop in rare cases where we keep getting the same index

while parent2_idx == parent1_idx
    
    while_loop_stop = while_loop_stop + 1;
    R = rand(); % in [0,1]
    
    if while_loop_stop > 20
        break;
    end
    
    for i = 1:nPop
        
        if R > cumsum(i)
            
            parent2_idx = i - 1;
            break;
        end
    end
end

p1Gene = tempGene(parent1_idx,:);
p2Gene = tempGene(parent2_idx,:);

end