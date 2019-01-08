function [wingDefs,foilData] = define_wing(options)
%% Define wing

Bezier = options.Bezier;
n = options.WingPartitions;

%% 2D Aerofoil Section Defintion

% Define sections as Bezier curves
if Bezier
    
    BezierControlPoints = options.BezierControlPoints;
    
    % Control point min/max coordinates
    minSec = [1, 0.7, 0.5,  0.3,    0.1,    0,      0;  % xu
        1,  0.7,    0.5,    0.3,    0.1,    0,      0;  % xl
        0,  0.015,  0.02,   0.05,   0.02,   0.05    0;  % zu
        0, -0.035, -0.04,  -0.07,  -0.04,  -0.05,   0]; % zl
    
    maxSec = [1, 0.9, 0.7,  0.5,    0.3,    0,      0;  % xu
        1,  0.9,    0.7,    0.5,    0.3,    0,      0;  % xl
        0,  0.035,  0.04,   0.07,   0.04,   0.05,   0;  % zu
        0, -0.015, -0.02,  -0.05,  -0.02,  -0.05,   0]; % zl
    
    [~,minControlPoints] = size(minSec);
    [~,maxControlPoints] = size(maxSec);
    
    if any([minControlPoints,maxControlPoints] ~=  BezierControlPoints)
        
        error("Number of Bezier control points not consistent with options definition")
    end
    
    % Transform matrix to single array
    minSec = reshape(minSec',1,[]);
    maxSec = reshape(maxSec',1,[]);
    
    minSec = repmat(minSec,1,n+1);
    maxSec = repmat(maxSec,1,n+1);
    
    foilData = length(minSec);
    
    sectionStr = "Bezier";
    sectionCons = "~";
    sectionTrans = "~";
else
    
    % Define sections be pre-loaded data files
    minSec = 1;
    
    % Load coordinates of 2D aerofoil sections into matrices within cell,
    % max defined by number of stored data files
    [foilData,maxSec] = getaerofoilsecdata();
    
    sectionStr = "Section";
    sectionCons = "Floor";
    sectionTrans = "~";
end

%   "Variables",            "VarMin",           "VarMax",   "Conditions"                "Transformations"
wingDefs = {...
    "Dihedral",             0,                  20,         "~",                        "~";
    "Chord",               [NaN,NaN,NaN,NaN],  [NaN,NaN,NaN,NaN],   "< Previous",               ".*AftLength";
    "LESweep",             [0,0,0],            [80,80,80],  "~"                         "~";
    "Semispan",            [0,0,0],            [5,5,5],    ["Minimum 0.5" "sum > 2"],   "~";
    sectionStr,             minSec,             maxSec,     sectionCons,                sectionTrans;
    "xOffset",             -0.25,               0.5,        "~",                        ".*AftLength";
    "zOffset",             -0.5,                0,          "~",                        ".*AftHalfHeight"};

[row,~] = size(wingDefs);

%% Error Checking
for i = 1:row
    
    name = wingDefs{i,1};
    
    if name == "Chord"
        
        if length(wingDefs{i,2}) ~= n+1 || length(wingDefs{i,3}) ~= n+1
            
            error("Chord min/max must be equal to number of partitions + 1 (n + 1)")
        end
    end
    
    if name == "LESweep"
        
        if length(wingDefs{i,2}) ~= n || length(wingDefs{i,3}) ~= n
            
            error("Leading edge sweep min/max must be equal to number of partitions (n)")
        end
    end
    
    if name == "Semispan"
        
        if length(wingDefs{i,2}) ~= n || length(wingDefs{i,3}) ~= n
            
            error("Semispan min/max must be equal to number of partitions (n)")
        end
    end
end