function [transformedPos] = hardtransform(initialPos,cond,varArray,direction)
%% Transformation from particle properties to physical or reverse
% eg. physicalPos(wingroot) = parPos(wingroot)*aftfuselagelength
% Direction can be set to "inverse" to reverse process
% Linear multiplications/division assumed (ratios)

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
    
    if exist('direction','var') && direction == "inverse"
        
        if contains(equation,"/")
            transformedPos(:,array) = target.*dimensionaliser;
        elseif contains(equation,"*")
            transformedPos(:,array) = target./dimensionaliser;
        end
        
    else
        
        if contains(equation,"/")
            transformedPos(:,array) = target./dimensionaliser;
        elseif contains(equation,"*")
            transformedPos(:,array) = target.*dimensionaliser;
        end
    end
end

end