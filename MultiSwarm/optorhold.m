function definition = optorhold(definition)

[row,col] = size(definition);

% Grab titles
for i=col:-1:1
    
    Headers(i) = definition{1,i};
end

defCol = Headers == "Variables";
minCol = Headers == "VarMin";
maxCol = Headers == "VarMax";
optCol = Headers == "Optimise/Hold";

% Now delete headers row
definition(1,:) = [];

for i = 1:row-1
    
    name = definition{i,defCol};
    optorhold = definition{i,optCol};
    varMin = definition{i,minCol};
    varMax = definition{i,maxCol};
    
    if optorhold == "Hold"
        
        if varMin == varMax
            
            fprintf("%s set at definition values as varMin = varMax \n", name)
        else

            definition{i,minCol}(:) = NaN;
            definition{i,maxCol}(:) = NaN;
            
            fprintf("%s will be set at baseline/standard constant values \n", name)
        end
    end
end
    