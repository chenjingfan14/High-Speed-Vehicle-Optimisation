function [condCell,varArray,var1,var2,nVar] = translate(variCons)
%% Conditioning initialisation
% Creates cell to be read for conditioning and indexing purposes

[dim1,dim2] = size(variCons);

preVar = 0;

condCell = cell(dim1-1,dim2);

var1 = [];
var2 = [];

foilData = getaerofoilsecdata();

[nFoils,~] = size(foilData);

% Grab titles
for i=dim2:-1:1
    Headers(i) = variCons{1,i};
end

% Should always exist regardless
nameCol = Headers == "Variables";

% One of should exist
numCol = Headers == "Num Of";

if any(numCol)
    getDim = false;
else
    numCol = Headers == "Values";
    getDim = true;
end

% May only exist for optimisation variables (not baseline or standard)
minCol = Headers == "VarMin";
maxCol = Headers == "VarMax";
conCol = Headers == "Conditions";
transCol = Headers == "Transformations";

%% VARMIN FOR OPT == VALUES FOR BASE/STANDARD. DO SOMETHING ABOUT IT

for i=2:dim1
    
    row = i - 1;
    
    name = variCons{i,nameCol};
    num = variCons{i,numCol};
        
    if isnumeric(num)
        
    elseif num == "~"
        
        num = 1;
        
    else
        
        for j = length(num):-1:1
            
            
            for k = 1:nFoils
                
                str = foilData{k,1};
                
                if contains(str,num(j))
                    newVal(j) = k;
                    break
                end
                
            end
        end
        
        num = newVal;
        
    end
    
    if getDim
        var = numel(num);
    else
        var = num;
    end
    
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
    if any(minCol) && any(maxCol)
        min = variCons{i,minCol};
        max = variCons{i,maxCol};
        
        if ~isnumeric(min)
            min = NaN;
        end
        
        if ~isnumeric(max)
            max = NaN;
        end
        
        min = repmat(min,1,var);
        max = repmat(max,1,var);
        
        var1 = [var1, min];
        var2 = [var2, max];
        
    else % varMin covers as value output
        
        var1 = [var1, num];
        
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
        if trans == "~"
            transFun = [];
        else
            transFun = trans;
        end
        
        condCell{row,6} = transFun;
        
    end
    
    preVar = nVar;
    
end

end
