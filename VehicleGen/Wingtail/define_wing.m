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
        0, -0.035, -0.04,  -0.07,  -0.04,  -0.05,   0]'; % zl
    
    maxSec = [1, 0.9, 0.7,  0.5,    0.3,    0,      0;  % xu
        1,  0.9,    0.7,    0.5,    0.3,    0,      0;  % xl
        0,  0.035,  0.04,   0.07,   0.04,   0.05,   0;  % zu
        0, -0.015, -0.02,  -0.05,  -0.02,  -0.05,   0]'; % zl
    
    % Use to hold sections at base/standard throughout simulation
%     [minSec,maxSec] = deal(nan(BezierControlPoints + 1,4));
    
    [minControlPoints,~] = size(minSec(2:end,:));
    [maxControlPoints,~] = size(maxSec(2:end,:));
    
    if any([minControlPoints,maxControlPoints] ~=  BezierControlPoints)
        
        error("Number of Bezier control points not consistent with options definition")
    end
    
    minSec = repmat(minSec, 1, 1, n + 1);
    maxSec = repmat(maxSec, 1, 1, n + 1);
    
    foilData = BezierControlPoints + 1;
    
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

%%
wingDefs = {...
    "Variables",            "VarMin",           "VarMax",   "Conditions",               "Transformations",      "Optimise/Hold"
    "Dihedral",             0,                  20,         "~",                        "~",                    "Optimise";
    "Chord",               [0.1,0.1,0.1,0.1],  [1,1,1,1],   "< Previous",               ".*AftLength",          "Optimise";
    "TESweep",             [-20,-20,-20],      [20,20,20],  "~"                         "~",                    "Optimise";
    "Semispan",            [0,0,0],            [5,5,5],    ["Minimum 0.5" "sum > 2"],   "~",                    "Optimise";
    sectionStr,             minSec,             maxSec,     sectionCons,                sectionTrans,           "Optimise";
    "xOffset",             -0.25,               0.5,        "~",                        ".*AftLength",          "Optimise";
    "zOffset",             -0.5,                0,          "~",                        ".*AftHalfHeight",      "Optimise"};

wingDefs = optorhold(wingDefs);

[row,~] = size(wingDefs);

%% Error Checking
for i = 1:row
    
    name = wingDefs{i,1};
    varMin = wingDefs{i,2};
    varMax = wingDefs{i,3};
    
    if name == "Chord"
        
        if length(varMin) ~= n+1 || length(varMax) ~= n+1
            
            error("Chord min/max must be equal to number of partitions + 1 (n + 1)")
        end
    end
    
    if name == "LESweep"
        
        if length(varMin) ~= n || length(varMax) ~= n
            
            error("Leading edge sweep min/max must be equal to number of partitions (n)")
        end
    end
    
    if name == "Semispan"
        
        if length(varMin) ~= n || length(varMax) ~= n
            
            error("Semispan min/max must be equal to number of partitions (n)")
        end
    end
end