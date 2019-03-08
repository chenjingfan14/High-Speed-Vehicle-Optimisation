%Written by Seyedali Mirjalili 

function [c1Gene , c2Gene] = crossover_continious(p1Gene , p2Gene, Pc, nVar, varMin, varMax)

c1Gene = zeros(1,nVar);
c2Gene = zeros(1,nVar);

for k = 1 : nVar
    
    beta = rand();
    c1Gene(k) = beta .* p1Gene(k) + (1-beta)*p2Gene(k); 
    c2Gene(k) = (1-beta) .* p1Gene(k) + beta*p2Gene(k);
    
    if c1Gene(k) > varMax(k) 
        c1Gene(k)  =  varMax(k);
    end
    
    if c1Gene(k) < varMin(k)
        c1Gene(k) = varMin(k);
    end
    
    if c2Gene(k) > varMax(k) 
        c2Gene(k)  =  varMax(k);
    end
    
    if c2Gene(k) < varMin(k)
        c2Gene(k) = varMin(k);
    end
end

R1 = rand();

if R1 > Pc

    c1Gene = p1Gene;
end

R2 = rand();

if R2 > Pc

    c2Gene = p2Gene;
end

end