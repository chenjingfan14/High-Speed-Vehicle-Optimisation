function [parPos,physicalPos] = conditioning(parPos,cond,varArray)
%% Function introduces conditions and ensures they are met
% Conditions placed on particle positions, translated to physical positions
% by ratios eg. physicalPos(wingroot) = parPos(wingroot)*aftfuselagelength

[dim,~] = size(cond);

dimArray = 1:dim;

%% Find variables corresponding to AftBodyHeight
bodyHeightCon = any(varArray == ["zUpperRad","SideLength","zLowerRad"]',1);
bodyHeightInd = zeros(sum(bodyHeightCon),1);

j = 1;
for i = dimArray(bodyHeightCon)
    bodyHeightInd(j) = cond{i,2};
    j = j + 1;
end

AftHeight = sum(parPos(:,bodyHeightInd),2); % Used in sym equations

%% Conditioning Loop
for i = 1:dim
    
    [targArray,consArray,settoArray,equation] = deal(cond{i,2:5});
    
    if ~isempty(equation)
        [~,numCons] = size(consArray);
        switch numCons
            case 0
                target = parPos(:,targArray); % Used in sym equation
                funAns = double(subs(str2sym(equation)));
                parPos(:,targArray) = funAns;
                
            case 1
                target = parPos(:,targArray);
                
                if contains(equation,"constraint")
                    constraint = parPos(:,consArray); % Used in sym equations
                    setto = parPos(:,settoArray);
                else
                    setto = settoArray;
                end
                
                comp = size(target);
                if ~isequal(size(setto),comp)
                    setto = repmat(setto,comp);
                end
                
                if contains(equation,["<",">","=="])
                    funAns = logical(subs(str2sym(equation)));
                    target(~funAns) = setto(~funAns);
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
                    
                    comp = size(target);
                    if ~isequal(size(setto),comp)
                        setto = repmat(setto,comp);
                    end
                    
                    if contains(equation,["<",">","=="])
                        funAns = logical(subs(str2sym(equation)));
                        target(~funAns) = setto(~funAns);
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
        physicalPos(:,array) = double(str2sym(equation));
    end
    
end
