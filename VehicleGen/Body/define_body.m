function bodyDefs = define_body(aft,fore,nose)
%% Define body

bodyDefs = [];

if aft
    
%       "Variables",            "VarMin",           "VarMax",   "Conditions"                "Transformations"
    aftDefs = {...
        "UpperLength",          NaN,                NaN,        "~",                        "~";
        "yUpperRad",            NaN,                NaN,        "~",                        "~";
        "yBotRatio",            NaN,                NaN,        "~",                        "~";
        "zUpperRad",            NaN,                NaN,        "~",                        "~";
        "SideLength",           NaN,                NaN,        "~",                        "~";
        "zLowerRad",            NaN,                NaN,        "~",                        "~";
        "AftLength",            NaN,                NaN,        "~",                        "~"};
    
    bodyDefs = [bodyDefs; aftDefs];
    
end

if fore
    
%       "Variables",            "VarMin",           "VarMax",   "Conditions"                "Transformations"
    foreDefs = {...
        "ForeLength",           NaN,                NaN,        "~",                        "~"};
    
    bodyDefs = [bodyDefs; foreDefs];
end

if nose
    
%       "Variables",            "VarMin",           "VarMax",   "Conditions"                "Transformations"
    noseDefs = {...
        "NoseRad",              NaN,                NaN,        "~",                        "~";
        "NoseLength",           NaN,                NaN,        "~",                        ".*NoseRad";
        "zNoseOffset",          NaN,                NaN,        "~",                        ".*AftHalfHeight"};
    
    bodyDefs = [bodyDefs; noseDefs];
end