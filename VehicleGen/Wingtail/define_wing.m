function [wingDefs] = define_wing(options)
%% Define wing

n = options.WingPartitions;

%%
wingDefs = {...
    "Variables",            "VarMin",           "VarMax",   "Conditions",               "Transformations",      "Optimise/Hold"
    "Dihedral",             0,                  20,         "~",                        "~",                    "Optimise";
    "Chord",               [0.1,0.1,0.1],      [1,1,1],     "< Previous",               ".*AftLength",          "Optimise";
    "TESweep",             [-20,-20],          [20,20],     "~"                         "~",                    "Optimise";
    "Semispan",            [2,1],              [5,5],      ["Minimum 0.5" "sum > 2"],   "~",                    "Optimise";
    "xOffset",             -0.25,               0.5,        "~",                        ".*AftLength",          "Hold";
    "zOffset",             -0.5,                0,          "~",                        ".*AftHalfHeight",      "Hold"};

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