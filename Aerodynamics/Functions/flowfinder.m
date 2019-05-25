function [points,newProperties,bodyPart,partType,impactMethod,shadowMethod] = flowfinder(properties,options)
%% Aerodynamic prediction method mixer
% Defines which prediction method to be used for each part and
% impact/shadow flow

% Ensure no cells are empty, if so delete them
properties = properties(~cellfun('isempty',properties));

dim = numel(properties);

numParts = ones(dim,1);

viscous = options.Viscous;
quad = options.Quad;

for i=1:dim
    
    propPoints = properties{i}.Points;
    
    if iscell(propPoints)
        
        numParts(i) = length(propPoints);
    end
end

%% Create points structure
% Count backwards for preallocation purposes
sumParts = sum(numParts);
id = sumParts;

for i = dim:-1:1
    
    name = properties{i}.Name;
    propPoints = properties{i}.Points;
    dim2 = numParts(i);
    
    for j = dim2:-1:1
        % Calculate additional panel properties and assign data to separate
        % points structure outwith properties cell
        if iscell(propPoints)
            
            points(id) = normals(propPoints{j},viscous,quad);
        else
            points(id) = normals(propPoints,viscous,quad);
        end
        
        if name == "aerofoil l"
            
            points(id).norm = -points(id).norm;
            points(id).unitNorm = -points(id).unitNorm;
        end
        
        %% DELETE? Need same amount of properties as point structures?
        newProperties{id} = properties{i};
         % Empty points in properties cell as we now have separate points struct
        newProperties{id}.Points = [];
        
        id = id - 1;
        
    end
end

bodyPart = true(sumParts,1);
[impactMethod,shadowMethod,partType] = deal(zeros(sumParts,1));

partType = string(partType);

%% Prediction method matrices
% Numbers correspond to switch case in parent function (aeroprediction)

for i=1:sumParts
    
    partType(i) = newProperties{i}.Name;
    
    % Prediction method mixer: Change method via aeroprediction
    % Impact:   1 - Modified Newtonian
    %           2 - Modified Newtonian + Prandtl-Meyer
    %           3 - Oblique Shock + Prandtl-Meyer
    %           4 - Tangent Wedge/Cone
    % Shadow:   1 - Newtonian/Base Pressure
    %           2 - Prandtl-Meyer
    switch partType(i)
        case {"aerofoil","wing","tail","aerofoil l","aerofoil u"}
            
            bodyPart(i) = false;
            impactMethod(i) = 3;
            shadowMethod(i) = 2;
            
        case "aftbody"
            
            impactMethod(i) = 3;
            shadowMethod(i) = 2;
            
        case "forebody"
            
            impactMethod(i) = 3;
            shadowMethod(i) = 2;
            
        case "nose"
           
            impactMethod(i) = 3;
            shadowMethod(i) = 2;
            
        case "test"
            
            impactMethod(i) = 3;
            shadowMethod(i) = 2;     
            
        otherwise
            
            bodyPart(i) = false;
            impactMethod(i) = 3;
            shadowMethod(i) = 2;    
    end 
end

if isfield(options,'ImpactMethod')
    
    impactMethod = options.ImpactMethod;
    shadowMethod = options.ShadowMethod;
end

end