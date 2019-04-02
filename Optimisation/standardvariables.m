function variCons = standardvariables(variCons,options)
%% Standardised variables to be held constant during optimisation

aerofoilMethod = options.AerofoilMethod;
wing = options.Wing;
aft = options.Aft;
fore = options.Fore;
nose = options.Nose;
control = options.Control;
baseline = options.Baseline;

%% 2D Aerofoil Section Defintion

if aerofoilMethod == "BP3434"
    
elseif aerofoilMethod == "Bezier"
    
    % Define sections as Bezier curves
    BezierControlPoints = options.BezierControlPoints;
    
    % Control point min/max coordinates
    standSec = [1,  0.9,    0.7,    0.5,    0.3,    0,      0;
        1,  0.9,    0.7,    0.5,    0.3,    0,      0;
        0,  0.035,  0.04,   0.07,   0.04,   0.05,   0;
        0, -0.015, -0.02,  -0.05,  -0.02,  -0.05,   0]';
    
    [standControlPoints,~] = size(standSec(2:end,:));
    
    if standControlPoints ~=  BezierControlPoints
        
        error("Number of Bezier control points not consistent with options definition")
    end
    
    sectionStr = "Bezier";
    
elseif aerofoilMethod == "BezierTC"
    
    BezierControlPoints = options.BezierControlPoints;
    
    % Control point min/max coordinates
    standSec = [1, 0.7,   0.5,    0.3,    0.1,    0,      0;  % xc
        1,      0.7,    0.5,    0.3,    0.1,    0,      0;  % xt
       -0.03,  -0.1,   -0.1,   -0.05,  -0.05,   0,      0;  % zc
        0,     -0.05,  -0.05,  -0.05,  -0.05,   0.05,   0]'; % zt
    
    [standSecPoints,~] = size(standSec(2:end,:));
    
    if standSecPoints ~=  BezierControlPoints
        
        error("Number of Bezier control points not consistent with options definition")
    end
    
    sectionStr = "Bezier";
    
elseif aerofoilMethod == "Preloaded"
    
    % Define sections be pre-loaded data files
    standSec = 1;
    sectionStr = "Section";
    
else
    error('No standard variables set for aerofoil creation method %s', aerofoilMethod);
end

Definition = {...
    "Variables",        "Standard Value"};

wingDefs = {...
    "Dihedral",         5;
    "Chord",            0.8;
    "TESweep",          0;
    "Semispan",         2;
    sectionStr,         standSec;
    "xOffset",          0;
    "zOffset",          0};

aftDefs = {...
    "UpperLength",      0;
    "yUpperRad",        0.88;
    "yBotRatio",        0.1;
    "zUpperRad",        0.88;
    "SideLength",       0.759;
    "zLowerRad",        0.05;
    "AftLength",        11.8956};

% Cylindrical Aftbody
% standardAft = {...
%     "UpperLength",      0;
%     "yUpperRad",        1;
%     "yBotRatio",        1;
%     "zUpperRad",        1;
%     "SideLength",       0;
%     "zLowerRad",        1;
%     "AftLength",        8};

foreDefs = {...
    "ForeLength",       4.4238};

noseDefs = {...
    "NoseRad",          0.155;
    "NoseLength",       0.7193;
    "zNoseOffset",     -0.6};

controlDefs = {...
    "ControlSpan",     [0.4,0.7];
    "ControlChord",     0.7};

[rows,~] = size(variCons);

if wing
    Definition = [Definition; wingDefs];
end

if aft
    Definition = [Definition; aftDefs];
end

if fore
    Definition = [Definition; foreDefs];
end

if nose
    Definition = [Definition; noseDefs];
end

if control
    Definition = [Definition; controlDefs];
end

if baseline
    
    baseConfig = options.Base;
    baseVariables = baseConfig.Variables;
    baseVarArray = baseConfig.VarArray;
    
    name = strings(rows,1);
    
    for i = 2:rows
        
        name(i) = variCons{i,1};
        
        varMin = variCons{i,2};
        varMax = variCons{i,3};
        
        minCon = isnan(varMin);
        maxCon = isnan(varMax);
        
        dims = size(varMin);
        elements = 1:numel(varMin);
        
        minArray = elements(minCon);
        maxArray = elements(maxCon);
        
        if any(minCon | maxCon)
            
            if name(i) ~= name(i-1)
                
                print = false;
            end
            
            stand = Definition{i,2};
            base = baseVariables(baseVarArray == name(i));
            
            baseDims = size(base);
            
            %% VarMin Loop
            use = min(dims(end),sum(minCon(:)));
            minArray = reshape(minArray,[],use);
            
            for j = 1:use
                
                ID = minArray(:,j);
                
                % If yes, sections wish to be set as stardard. Must check that
                % optimisation variables and base variables are the same for the wing
                % (in terms of nPartitions) otherwise will not work
                if j <= baseDims(end)
                    
                    varMin(ID) = base(ID);
                else
                    
                    if ~print
                        
                        if j == 1
                            
                            fprintf('Standard %s used as base definition not applicable \n', name(i));
                        else
                            
                            fprintf('Base %s used outboard until partition %i (or section %i) \n', name(i), j-1, j)
                        end
                        
                        print = true;
                    end
                    
                    varMin(ID) = stand;
                end
            end
            
            variCons{i,2} = varMin;
            
            %% VarMax Loop
            use = min(dims(end),sum(maxCon(:)));
            maxArray = reshape(maxArray,[],use);
            
            for j = 1:use
                
                ID = maxArray(:,j);
                
                % If yes, sections wish to be set as stardard. Must check that
                % optimisation variables and base variables are the same for the wing
                % (in terms of nPartitions) otherwise will not work
                if j <= baseDims(end)
                    
                    varMax(ID) = base(ID);
                else
                    
                    if ~print
                        
                        if j == 1
                            
                            fprintf('Standard %s used as base definition not applicable \n', name(i));
                        else
                            
                            fprintf('Base %s used outboard until partition %i (or section %i) \n', name(i), j-1, j)
                        end
                        
                        print = true;
                    end
                    
                    varMax(ID) = stand;
                end
            end
            
            variCons{i,3} = varMax;
        end
    end
else
    
    % If no baseline specific, replace all desired variables with standard
    % ones specified above
    for i = 2:rows
        
        varMin = variCons{i,2};
        varMax = variCons{i,3};
        
        stand = Definition{i,2};
        
        minCon = isnan(varMin);
        maxCon = isnan(varMax);
        
        if length(stand) == 1
            varMin(minCon) = stand;
            varMax(maxCon) = stand;
        else
            varMin(minCon) = stand(minCon);
            varMax(maxCon) = stand(maxCon);
        end
        
        variCons{i,2} = varMin;
        variCons{i,3} = varMax;
    end
    
end

end