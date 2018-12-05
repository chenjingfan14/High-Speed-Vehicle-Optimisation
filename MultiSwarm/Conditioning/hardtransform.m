function [transformedPos] = hardtransform(initialPos,cond,varArray)
%% Transformation from particle properties to physical
% eg. physicalPos(wingroot) = parPos(wingroot)*aftfuselagelength

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
% From particle position to physical configuration inputs
for i = 1:dim
    
    [array,equation] = deal(cond{i,[2 6]});
    
    if ~isempty(equation)
        switch i
            
            % Transform
            % Reverse
            
            case 2
                transformedPos(:,array) = initialPos(:,array).*AftLength;
                % transformedPos(:,array) = initialPos(:,array)./AftLength;
            case 6
                transformedPos(:,array) = initialPos(:,array).*AftLength;
                % transformedPos(:,array) = initialPos(:,array)./AftLength;
            case 7
                transformedPos(:,array) = initialPos(:,array).*AftHalfHeight;
                % transformedPos(:,array) = initialPos(:,array)./AftHalfHeight;
            case 16
                transformedPos(:,array) = initialPos(:,array).*NoseRad;
                % transformedPos(:,array) = initialPos(:,array)./NoseRad;
            case 17
                transformedPos(:,array) = initialPos(:,array).*AftHalfHeight;
                % transformedPos(:,array) = initialPos(:,array)./AftHalfHeight;
        end
    end
    
end
