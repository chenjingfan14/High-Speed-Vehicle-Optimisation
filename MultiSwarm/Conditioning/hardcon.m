function [parPos,physicalPos] = hardcon(parPos,nPop,cond,varArray)
%% Function introduces conditions and ensures they are met
% Conditions placed on particle positions, translated to physical positions
% by ratios eg. physicalPos(wingroot) = parPos(wingroot)*aftfuselagelength

[dim,~] = size(cond);

dimArray = 1:dim;
popArray = (1:nPop)';
%% Find variables corresponding to AftBodyHeight, etc
bodyHeightCon = any(varArray == ["zUpperRad","SideLength","zLowerRad"]',1);
bodyHeightInd = zeros(sum(bodyHeightCon),1);

j = 1;
for i = dimArray(bodyHeightCon)
    bodyHeightInd(j) = cond{i,2};
    j = j + 1;
end

AftHeight = sum(parPos(:,bodyHeightInd),2)/2; % Used in sym equations

aftLengthCon = varArray == "AftLength";
aftLengthInd = cond{dimArray(aftLengthCon),2};
AftLength = parPos(:,aftLengthInd);

noseRadCon = varArray == "NoseRad";
noseRadInd = cond{dimArray(noseRadCon),2};
NoseRad = parPos(:,noseRadInd);

%% Conditioning Loop
for i = 1:dim
    
    [targArray,consArray,settoArray,equation] = deal(cond{i,2:5});
    
    ID = popArray + (targArray-1)*nPop;
    
    if ~isempty(equation)
        [numConditions,numConstraints] = size(consArray);
        
        for ii = 1:numConditions
            
            eq = equation(ii);
            
            switch numConstraints
                case 0
                    target = parPos(ID);
                    
                    if contains(eq,"floor")
                        funAns = floor(target);
                    elseif contains(eq,"ceil")
                        funAns = ceil(target);
                    end
                    
                    parPos(ID) = funAns;
                    
                case 1
                    target = parPos(ID);
                    
                    if contains(eq,"constraint")
                        constraint = parPos(:,consArray);
                        setto = parPos(:,settoArray);
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
                    
                    parPos(ID(~satisfied)) = setto(~satisfied);
                    
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
                        
                        target = parPos(:,j);
                        
                        if contains(eq,"constraint")
                            constraint = parPos(:,consArray(counter)); % Used in sym equations
                            setto = parPos(:,settoArray(counter));
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
                        
                        parPos(:,j) = target;
                        counter = counter + 1;
                        
                    end
            end
        end
    end
end

physicalPos = parPos;

%% Alterations
% From particle to position to physical configuration inputs
for i = 1:dim
    
    [array,equation] = deal(cond{i,[2 6]});
    
    if ~isempty(equation)
        switch i
            case 2
                physicalPos(:,array) = parPos(:,array).*AftLength;
            case 6
                physicalPos(:,array) = parPos(:,array).*AftLength;
            case 7
                physicalPos(:,array) = parPos(:,array).*AftHeight;
            case 16
                physicalPos(:,array) = parPos(:,array).*NoseRad;
            case 17
                physicalPos(:,array) = parPos(:,array).*AftHeight;
        end
    end
    
end
