function [points,numParts,bodyPart,impactMethod,shadowMethod] = flowfinder(properties)
%% Aerodynamic prediction method mixer
% Defines which prediction method to be used for each part and
% impact/shadow flow

% Ensure no cells are empty, if so delete them
properties = properties(~cellfun('isempty',properties));

dim = length(properties);

bodyPart = true(dim,1);
[numParts,impactMethod,shadowMethod] = deal(zeros(dim,1));

%% Prediction method matrices
% Numbers correspond to switch case in parent function (aeroprediction)
for i=1:dim
    propPoints = properties{i}.Points;
    numParts(i) = length(propPoints);
    type = propPoints.Name;
    
    % Prediction method mixer: Change method via aeroprediction
    % Impact:   1 - Modified Newtonian
    %           2 - Modified Newtonian + Prandtl-Meyer
    %           3 - Oblique Shock + Prandtl-Meyer
    %           4 - Tangent Wedge/Cone
    % Shadow:   1 - Newtonian/Base Pressure
    %           2 - Prandtl-Meyer
    switch type
        case "aerofoil"
            
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
        id = id - 1;
        
    end
    
    % Empty points in properties cell as we now have separate points struct
    properties{i}.Points = [];
end