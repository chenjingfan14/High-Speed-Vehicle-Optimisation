function [parPos,physicalPos] = hardcon(parPos,cond,varArray)
%% Function introduces conditions and ensures they are met
% Conditions placed on particle positions, translated to physical positions
% by ratios eg. physicalPos(wingroot) = parPos(wingroot)*aftfuselagelength

[dim,~] = size(cond);

dimArray = 1:dim;

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
    
    if ~isempty(equation)
        [~,numCons] = size(consArray);
        switch numCons
            case 0
                target = parPos(:,targArray);
                
                if contains(equation,"floor")
                    funAns = floor(target);
                elseif contains(equation,"ceil")
                    funAns = ceil(target);
                end
                
                parPos(:,targArray) = funAns;
                
            case 1
                target = parPos(:,targArray);
                
                if contains(equation,"constraint")
                    constraint = parPos(:,consArray); % Used in sym equations
                    setto = parPos(:,settoArray);
                else
                    constraint = consArray;
                    setto = settoArray;
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
                
                if contains(equation,"<")
                    satisfied = target < constraint;
                    target(~satisfied) = setto(~satisfied);
                elseif contains(equation,">")
                    satisfied = target > constraint;
                    target(~satisfied) = setto(~satisfied);
                end
                
                parPos(:,targArray) = target;
                
            otherwise
                counter = 1;
                for j = targArray
                    target = parPos(:,j);
                    
                    if contains(equation,"constraint")
                        constraint = parPos(:,consArray(counter)); % Used in sym equations
                        setto = parPos(:,settoArray(counter));
                    else
                        setto = settoArray(:,counter);
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
                    
                    if contains(equation,"<")
                        satisfied = target <= constraint;
                        target(~satisfied) = setto(~satisfied);
                    elseif contains(equation,">")
                        satisfied = target >= constraint;
                        target(~satisfied) = setto(~satisfied);
                    end
                    
                    parPos(:,j) = target;
                    counter = counter + 1;
                    
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