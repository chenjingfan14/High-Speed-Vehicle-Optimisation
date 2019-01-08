%% Define baseline body

bodyDefs = [];

if aft
    
    aftDefs = {...
        "UpperLength",  0;
        "yUpperRad",    0.88;
        "yBotRatio",    0.1;
        "zUpperRad",    0.88;
        "SideLength",   0.759;
        "zLowerRad",    0.05;
        "AftLength",    11.8956};
    
    bodyDefs = [bodyDefs; aftDefs];
end

if fore
    
    foreDefs = {...
        "ForeLength",   4.4238};
    
    bodyDefs = [bodyDefs; foreDefs];
end

if nose
    
    noseDefs = {...
        "NoseRad",      0.155;
        "NoseLength",   0.1115;
        "zNoseOffset", -0.6};

%     noseDefs = {...
%         "NoseRad",      0.155;
%         "NoseLength",   0.7194;
%         "zNoseOffset", -0.7105};
    
    bodyDefs = [bodyDefs; noseDefs];
end