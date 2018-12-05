function [transformedPos] = versatiletransform(initialPos,cond,varArray)
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

AftHalfHeight = sum(initialPos(:,bodyHeightInd),2)/2; % Used in sym equations

aftLengthCon = varArray == "AftLength";
aftLengthInd = cond{dimArray(aftLengthCon),2};
AftLength = initialPos(:,aftLengthInd);

noseRadCon = varArray == "NoseRad";
noseRadInd = cond{dimArray(noseRadCon),2};
NoseRad = initialPos(:,noseRadInd);

transformedPos = initialPos;
%% Alterations
% From particle to position to physical configuration inputs
for i = 1:dim
    
    [array,equation] = deal(cond{i,[2 6]});
    
    if ~isempty(equation)
        transformedPos(:,array) = double(str2sym(equation));
    end
    
end
