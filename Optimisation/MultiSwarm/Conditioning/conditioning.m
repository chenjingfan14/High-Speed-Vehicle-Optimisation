function [initialPos] = conditioning(initialPos,nPop,cond)
%% Function introduces conditions and ensures they are met
% Conditions placed on particle positions

[dim,~] = size(cond);

popArray = (1:nPop)';


%% Conditioning Loop
for i = 1:dim
    
    [targArray,consArray,settoArray,equation] = deal(cond{i,2:5});
    
    ID = popArray + (targArray-1)*nPop;
    
    if ~isempty(equation)
        
        [numConditions,numConstraints] = size(consArray);
        
        if all(isnan(consArray),2)
            numConstraints = 0;
        end
        
        for ii = 1:numConditions
            
            eq = equation(ii);
            
            switch numConstraints
                case 0
                    target = initialPos(ID);
                    
                    if contains(eq,"floor")
                        funAns = floor(target);
                    elseif contains(eq,"ceil")
                        funAns = ceil(target);
                    end
                    
                    initialPos(ID) = funAns;
                    
                case 1
                    target = initialPos(ID);
                    
                    if contains(eq,"constraint")
                        constraint = initialPos(:,consArray);
                        setto = initialPos(:,settoArray);
                    else
                        constraint = consArray(ii);
                        setto = settoArray(ii);
                    end
                    
                    % Sums all target variables, returns only 1 value,
                    % therefore targArray set to only 1 column (first)
                    if contains(eq,"sum")
                        target = sum(target,2);
                        ID = ID(:,1);
                    end
                    
                    [rows,cols] = size(target);
                    if ~isequal(size(setto,1),rows)
                        setto = repmat(setto,rows,1);
                    end
                    if ~isequal(size(constraint,1),rows)
                        constraint = repmat(constraint,rows,1);
                    end
                    
                    if ~isequal(size(setto,2),cols)
                        setto = repmat(setto,1,cols);
                    end
                    if ~isequal(size(constraint,2),cols)
                        constraint = repmat(constraint,1,cols);
                    end
                    
                    if contains(eq,"<")
                        satisfied = target < constraint;
                    elseif contains(eq,">")
                        satisfied = target > constraint;
                    end
                    
                    initialPos(ID(~satisfied)) = setto(~satisfied);
                    
                otherwise
                    counter = 1;
                    
                    % For occasions where multiple conditions are present.
                    % Currently these occasions are assumed only to be when
                    % such conditions are based on other particle parameters,
                    % ie. can only be indexes, not arbitrary number constraints
                    for j = targArray
                        
                        if isnan(consArray(counter))
                            counter = counter + 1;
                            continue
                        end
                        
                        target = initialPos(:,j);
                        
                        if contains(eq,"constraint")
                            constraint = initialPos(:,consArray(counter)); % Used in sym equations
                            setto = initialPos(:,settoArray(counter));
                        else
                            setto = settoArray(:,counter);
                        end
                        
                        if contains(eq,"sum")
                            target = sum(target,2);
                            targArray = targArray(1);
                        end
                        
                        [rows,cols] = size(target);
                        if ~isequal(size(setto,1),rows)
                            setto = repmat(setto,rows,1);
                        end
                        if ~isequal(size(constraint,1),rows)
                            constraint = repmat(constraint,rows,1);
                        end
                        
                        if ~isequal(size(setto,2),cols)
                            setto = repmat(setto,1,cols);
                        end
                        if ~isequal(size(constraint,2),cols)
                            constraint = repmat(constraint,1,cols);
                        end
                        
                        if contains(eq,"<")
                            satisfied = target <= constraint;
                            target(~satisfied) = setto(~satisfied);
                        elseif contains(eq,">")
                            satisfied = target >= constraint;
                            target(~satisfied) = setto(~satisfied);
                        end
                        
                        initialPos(:,j) = target;
                        counter = counter + 1;
                        
                    end
            end
        end
    end
end
