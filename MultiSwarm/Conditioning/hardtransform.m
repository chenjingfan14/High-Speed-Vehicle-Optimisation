function [transformedPos] = hardtransform(initialPos,cond,varArray)
%% Transformation from particle properties to physical
% eg. physicalPos(wingroot) = parPos(wingroot)*aftfuselagelength

[dim,~] = size(cond);

%% Find variables corresponding to AftBodyHeight, etc
bodyHeightCon = any(varArray == ["zUpperRad","SideLength","zLowerRad"],2);
AftHeight = sum(initialPos(:,bodyHeightCon),2);

aftLengthCon = varArray == "AftLength";
AftLength = initialPos(:,aftLengthCon);

noseRadCon = varArray == "NoseRad";
NoseRad = initialPos(:,noseRadCon);

transformedPos = initialPos;

%% Alterations
% From particle position to physical configuration inputs
for i = 1:dim
    
    [array,equation] = deal(cond{i,[2 6]});
    
    if isempty(equation)
        continue
    end
    
    target = initialPos(:,array);
    
    if contains(equation,"AftLength")
        dimensionaliser = AftLength;
    elseif contains(equation,"AftHeight")
        dimensionaliser = AftHeight;
    elseif contains(equation,"AftHalfHeight")
        dimensionaliser = AftHeight/2;
    elseif contains(equation,"NoseRad")
        dimensionaliser = NoseRad;
    end
    
    if contains(equation,"/")
        transformedPos(:,array) = target./dimensionaliser;
    elseif contains(equation,"*") 
        transformedPos(:,array) = target.*dimensionaliser;
    end
    
end

end