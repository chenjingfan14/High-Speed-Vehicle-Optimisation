function bodyDefs = define_body(aft,fore,nose)
%% Define body

bodyDefs = [];

if aft
    
    aftDefs = {...
        "Variables",        "VarMin",   "VarMax",   "Conditions",	"Transformations",  "Optimise/Hold"
        "UpperLength",      NaN,        NaN,        "~",            "~",                "Hold";
        "yUpperRad",        NaN,        NaN,        "~",            "~",                "Hold";
        "yBotRatio",        NaN,        NaN,        "~",            "~",                "Hold";
        "zUpperRad",        NaN,        NaN,        "~",            "~",                "Hold";
        "SideLength",       NaN,        NaN,        "~",            "~",                "Hold";
        "zLowerRad",        NaN,        NaN,        "~",            "~",                "Hold";
        "AftLength",        NaN,        NaN,        "~",            "~",                "Hold"};
    
    bodyDefs = [bodyDefs; aftDefs];
    
end

if fore
    
    foreDefs = {...
        "ForeLength",       NaN,        NaN,        "~",            "~",                "Hold"};
    
    bodyDefs = [bodyDefs; foreDefs];
end

if nose
    
    noseDefs = {...
        "NoseRad",          NaN,        NaN,        "~",            "~",                "Hold";
        "NoseLength",       NaN,        NaN,        "~",            ".*NoseRad",        "Hold";
        "zNoseOffset",      NaN,        NaN,        "~",            ".*AftHalfHeight",  "Hold"};
    
    bodyDefs = [bodyDefs; noseDefs];
end

bodyDefs = optorhold(bodyDefs);