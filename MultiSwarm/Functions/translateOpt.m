function [condCell,varArray,var1,var2,nVar] = translateOpt(variCons)
%% Conditioning initialisation
% Creates cell to be read for conditioning and indexing purposes

[rows,cols] = size(variCons);
colArray = 1:cols;

preVar = 0;

condCell = cell(rows-1,6);

var1 = [];
var2 = [];

% Grab titles
for i=cols:-1:1
    
    Headers(i) = variCons{1,i};
end

% Should always exist regardless
nameCol = Headers == "Variables";
valCol = any(Headers == ["VarMin","VarMax","Values"]');

valColArray = colArray(valCol);

% May only exist for optimisation variables (not baseline or standard)
conCol = Headers == "Conditions";
transCol = Headers == "Transformations";

for i=2:rows
    
    row = i - 1;
    
    name = variCons{i,nameCol};
    val = variCons{i,valCol};
    
    var = length(val);
    
    nVar = preVar + var;
    target = (preVar + 1) : nVar;
    
    if preVar == 0
        varArray(1,:) = name;
    else
        varArray(target,:) = name;
    end
    
    condCell{row,1} = name;
    condCell{row,2} = target;
    
    %% Min/Max Values
    if sum(valCol) > 1
        
        min = variCons{i,valColArray(1)};
        max = variCons{i,valColArray(2)};
        
        var1 = [var1, min];
        var2 = [var2, max];
    else % varMin covers as value output
        
        var1 = [var1, val];
    end
    
    %% Condition
    if any(conCol)
        
        [Constraint,Setto,conFun] = deal([]);
        con = variCons{i,conCol};
        
        if con ~= "~"
            
            for ii = numel(con):-1:1
                
                coni = con(ii);
                
                if contains(coni,"<")
                    condition = ' < ';
                end
                if contains(coni,">")
                    condition = ' > ';
                end
                if contains(coni,"==")
                    condition = ' == ';
                end
                if contains(coni,"sum")
                    constraint = regexp(coni,'[+-]?\d+\.?\d*', 'match');
                    setto = constraint;
                    equation = ['sum(target)' condition char(constraint)];
                    constraint = double(constraint);
                end
                if contains(coni,"Previous")
                    constraint = [NaN, target(1:end-1)];
                    setto = constraint;
                    equation = ['target' condition 'constraint'];
                end
                if contains(coni,"Next")
                    constraint = [target(2:end), NaN];
                    setto = constraint;
                    equation = ['target' condition 'constraint'];
                end
                if contains(coni,"Minimum")
                    
                    constraint = regexp(coni,'[+-]?\d+\.?\d*', 'match');
                    setto = 0;
                    equation = ['target > ' char(constraint)];
                    constraint = double(constraint);
                end
                if contains(coni,"Maximum")
                    
                    constraint = regexp(coni,'[+-]?\d+\.?\d*', 'match');
                    equation = ['target < ' char(constraint)];
                    constraint = double(constraint);
                end
                if contains(coni,"Floor")
                    equation = 'floor(target)';
                    constraint = NaN;
                    setto = NaN;
                end
                
                Constraint(ii,:) = constraint;
                Setto(ii,:) = setto;
                conFun{ii} = equation;
                
            end
            
        end
        
        condCell{row,3} = Constraint;
        condCell{row,4} = Setto;
        condCell{row,5} = string(conFun);
        
    end
    
    %% Transformation
    if any(transCol)
        
        trans = variCons{i,transCol};
        
        if trans ~= "~"
            
            transFun = trans;
            condCell{row,6} = transFun;
        end
    end
    
    preVar = nVar;
    
end

if numel(varArray) ~= numel(var1)
    error("Number of variables does not match corresponding variable names")
end

end
