function [points,newProperties,bodyPart,partType,impactMethod,shadowMethod] = flowfinder(properties)
%% Aerodynamic prediction method mixer
% Defines which prediction method to be used for each part and
% impact/shadow flow

% Ensure no cells are empty, if so delete them
properties = properties(~cellfun('isempty',properties));

dim = numel(properties);

numParts = zeros(dim,1);

for i=1:dim
    numParts(i) = length(properties{i}.Points);
end

%% Create points structure
% Count backwards for preallocation purposes
id = sum(numParts);

for i=dim:-1:1
    
    propPoints = properties{i}.Points;
    dim2 = length(propPoints);
    
    for j = dim2:-1:1
        % Calculate additional panel properties and assign data to separate
        % points structure outwith properties cell
        points(id) = normals(propPoints(j));
        
        newProperties{id} = properties{i};
%         newProperties{id}.Name = points.Name;
         % Empty points in properties cell as we now have separate points struct
        newProperties{id}.Points = [];
        
        id = id - 1;
        
    end
end

dim = sum(numParts);

bodyPart = true(dim,1);
[impactMethod,shadowMethod,partType] = deal(zeros(dim,1));

partType = string(partType);

%% Prediction method matrices
% Numbers correspond to switch case in parent function (aeroprediction)
for i=1:dim
    partType(i) = newProperties{i}.Name;
%     partType(i) = propPoints.Name;
    
    % Prediction method mixer: Change method via aeroprediction
    % Impact:   1 - Modified Newtonian
    %           2 - Modified Newtonian + Prandtl-Meyer
    %           3 - Oblique Shock + Prandtl-Meyer
    %           4 - Tangent Wedge/Cone
    % Shadow:   1 - Newtonian/Base Pressure
    %           2 - Prandtl-Meyer
    switch partType(i)
        case {"aerofoil","wing","tail"}
            
            bodyPart(i) = false;
            impactMethod(i) = 1;
            shadowMethod(i) = 1;
            
        case "aftbody"
            
            impactMethod(i) = 1;
            shadowMethod(i) = 1;
            
        case "forebody"
            
            impactMethod(i) = 1;
            shadowMethod(i) = 1;
            
        case "nose"
           
            impactMethod(i) = 1;
            shadowMethod(i) = 1;
            
        case "test"
            
            impactMethod(i) = 1;
            shadowMethod(i) = 1;
            
    end 
end

end