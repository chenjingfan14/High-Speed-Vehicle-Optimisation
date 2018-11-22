function [condCell,varArray,nVar] = translate(variCons)
%% Conditioning initialisation
% Creates cell to be read for conditioning and indexing purposes

[dim1,dim2] = size(variCons);

for i=dim1:-1:2
    varArray(i-1) = variCons{i,1};
    numArray(i-1) = i;
end

condCell = cell(dim1-1,dim2+2);
strcnt = 1;

for i=2:dim1
    
    name = variCons{i,1};
    num = variCons{i,2};
    con = variCons{i,3};
    
    if isnumeric(num)
        a = num*length(name);
    else
        a = length(name);
    end
    
    target = strcnt:strcnt+a-1;
    [Constraint,Setto,conFun] = deal([]);
    
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
            end
            
            Constraint(ii,:) = constraint;
            Setto(ii,:) = setto;
            conFun{ii} = equation;
            
        end
        
    end
    j = i - 1;
    
    condCell{j,1} = name;
    condCell{j,2} = target;
    condCell{j,3} = Constraint;
    condCell{j,4} = Setto;
    condCell{j,5} = string(conFun);
    strcnt = strcnt + a;
end

nVar = strcnt - 1;

for i = j:-1:1
    
    num = condCell{i,2};
    relation = variCons{i+1,4};
    condCell{i,6} = [];
    
    if relation ~= "~"
        
        if length(num) > 1
            parPos = ['parPos(:,[' num2str(num) '])'];
        else
            parPos = ['parPos(:,' num2str(num) ')'];
        end
        
        equation = ['(' parPos relation ')'];
        vars = allwords(relation);
        
        for ii = length(vars):-1:1
            
            var = char(vars(ii));
            row = contains(varArray,var);
            
            if any(row)
                
                tranNum = numArray(row)-1;
                tranNum = condCell{tranNum,2};
                str = ['physicalPos(:,' num2str(tranNum) ')'];
                equation = char(strrep(equation,var,str));
                
            else
                equation = char(strrep(equation,var,[var '(:,1)']));
            end
            
        end
        
        condCell{i,6} = equation;
        
    end
    
end
