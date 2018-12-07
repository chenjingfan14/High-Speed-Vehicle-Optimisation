function [condCell,vals,varArray] = translateVals(variCons)
%% Conditioning initialisation
% Creates cell to be read for conditioning and indexing purposes

[dim1,~] = size(variCons);

preVar = 0;
vals = [];
for i=2:dim1
    
    name = variCons{i,1};
    val = variCons{i,2};
    
    
    
    if isnumeric(val)
        
    else
        
        for j = length(val):-1:1
           
            
            for k = 1:nFoils
                
                str = foilData{k,1};
                
                if contains(str,val(j))
                    newVal(j) = k;
                    break
                end

            end
        end
        
        val = newVal;
        
    end
    
    vals = [vals, val];
    
    var = numel(val);
    
    nVar = preVar + var;
    target = (preVar + 1) : nVar;
    
    if preVar == 0
        varArray(1,:) = name;
    else
        varArray(target,:) = name;
    end
    
    preVar = nVar;
    
end

end
