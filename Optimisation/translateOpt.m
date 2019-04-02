function [condCell,varArray,var1,var2,nVar] = translateOpt(variCons)
%% Conditioning initialisation
% Creates cell to be read for conditioning and indexing purposes

[rows,cols] = size(variCons);

preVar = 0;

condCell = cell(rows-1,6);

var1 = 0;
var2 = 0;

% Grab titles
for i=cols:-1:1
    
    Headers(i) = variCons{1,i};
end

% Should always exist regardless
nameCol = Headers == "Variables";
valCol = any(Headers == ["VarMin","VarMax","Values"]');

% May only exist for optimisation variables (not baseline or standard)
conCol = Headers == "Conditions";
transCol = Headers == "Transformations";

for i = 2:rows
    
    row = i - 1;
    
    name = variCons{i,nameCol};
    
    valCell = variCons(i,valCol);
    
    %% Transform indidvidual characteristic values into 1D all encompassing arrays
    
    % Counter can only equal 1 or 2
    
    % If optimisation variCons being translated: 
    % var1 = varMin, var2 = varMax
    
    % If baseline variCons being translated:
    % var1 = variables, var2 not used
    for j = 1:numel(valCell)
        
        val = valCell{j};
        dims = size(val);
        
        % Need 1D arrays, so transform any higher dimension matrices into
        % such
        if numel(dims) == 3 && all(dims > 1) % Turns 3D into 1D
            
            val = reshape(val,1,[],dims(3));
            val = reshape(val,1,[]);
            
        elseif numel(dims) == 2 && all(dims > 1) % Turns 2D into 1D
            
            val = reshape(val,1,[]);
        end
        
        % Assign names to variables and form min and max arrays
        if j == 1
            
            % Find where current variables lie in complete variable array
            var = length(val);
            nVar = preVar + var;
            target = (preVar + 1) : nVar;
            
            % Assign current variables to complete variable array in
            % correct position    
            varArray(target,:) = name;
            var1(target,:) = val;
  
        elseif j == 2
            
            % As above however target has already been calculated, and
            % varArray has been assigned already. Here only assigning
            % complete maximum variable array (if applicable)
            var2(target,:) = val;
        end
    end
    
    condCell{row,1} = name;
    condCell{row,2} = target;
    
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
