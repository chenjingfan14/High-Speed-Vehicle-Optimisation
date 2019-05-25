function [foilDefs,options] = define_aerofoil(options)

aerofoilMethod = options.AerofoilMethod;
n = options.WingPartitions;

%% 2D Aerofoil Section Defintion

% Define sections as Bezier curves
switch aerofoilMethod
    
    case {"Parsec","PARSEC"}
        
        foilTitles = {...
            "Variables",    "VarMin",   "VarMax",   "Conditions",	"Transformations",	"Optimise/Hold"};
        
        foilDefs = {...
            "rleu",         0.01,       0.1,        "~",            "~",                "Optimise";          
            "rlel",        -0.1,       -0.01,       "~",            "~",                "Optimise";
            "xup",          0.2,        0.8,        "~",            "~",                "Optimise";
            "zup",          0.05,       0.15,       "~",            "~",                "Optimise";
            "zxxup",       -0.6,       -0.4,        "~",            "~",                "Optimise";
            "xlo",          0.2,        0.8,        "~",            "~",                "Optimise";
            "zlo",         -0.15,      -0.05,       "~",            "~",                "Optimise";
            "zxxlo",        0.4,        0.6,        "~",            "~",                "Optimise";
            "zte",         -0.05,       0.05,       "~",            "~",                "Optimise";
            "dzte",         0,          0.05,       "~",            "~",                "Optimise";
            "ate",         -10,         10,         "~",            "~",                "Optimise";
            "bte",          0,          20,         "~",            "~",                "Optimise"};
            
        foilDefs = [foilTitles; repmat(foilDefs, n + 1, 1)];
        
    case "BP3434"
        
        foilDefs = {...
            "Variables",    "VarMin",   "VarMax",   "Conditions",	"Transformations",	"Optimise/Hold";
            "gle",          0.05,       0.1,        "~",            "~",                "Optimise";
            "b0",           0.01,       0.1,        "~",            "~",                "Optimise";
            "b2",           0.1,        0.3,        "~",            "~",                "Optimise";
            "xc",           0.2,        0.5,        "~",            "~",                "Optimise";
            "yc",           0,          0.2,        "~",            "~",                "Optimise";
            "b17",          0,          0.9,        "~",            "~",                "Optimise";
            "zte",          0,          0.01,       "~",            "~",                "Optimise";
            "ate",          0.05,       0.1,        "~",            "~",                "Optimise";
            "rle",         -0.04,      -0.001,      "~",            "~",                "Optimise";
            "b8",           0,          0.7,        "~",            "~",                "Optimise";
            "xt",           0.15,       0.4,        "~",            "~",                "Optimise";
            "yt",           0.05,       0.15,       "~",            "~",                "Optimise";
            "b15",          0,          0.9,        "~",            "~",                "Optimise";
            "dzte",         0,          0.001,      "~",            "~",                "Optimise";
            "bte",          0.001,      0.3,        "~",            "~",                "Optimise";
            "b1",           0.08,       0.17,       "~",            "~",                "Optimise"};
        
        [row,~] = size(foilDefs);
        
        for i = 2:row
            
            foilDefs{i,2} = repmat(foilDefs{i,2}, 1, n + 1);
            foilDefs{i,3} = repmat(foilDefs{i,3}, 1, n + 1);
        end
        
    case "Bezier"
        
        BezierControlPoints = options.BezierControlPoints;
        
        % Control point min/max coordinates
        %     minSec = [1, 0.7, 0.5,  0.3,    0.1,    0,      0;  % xu
        %         1,  0.7,    0.5,    0.3,    0.1,    0,      0;  % xl
        %         0,  0.015,  0.02,   0.05,   0.02,   0.05    0;  % zu
        %         0, -0.035, -0.04,  -0.07,  -0.04,  -0.05,   0]'; % zl
        %
        %     maxSec = [1, 0.9, 0.7,  0.5,    0.3,    0,      0;  % xu
        %         1,  0.9,    0.7,    0.5,    0.3,    0,      0;  % xl
        %         0,  0.035,  0.04,   0.07,   0.04,   0.05,   0;  % zu
        %         0, -0.015, -0.02,  -0.05,  -0.02,  -0.05,   0]'; % zl
        
        minSec = [1, 0.7, 0.5,  0.3,    0.1,    0,      0;  % xu
            1,      0.7,    0.5,    0.3,    0.1,    0,      0;  % xl
            0,      0,      0,      0,      0,      0       0;  % zu
            -0.03,  -0.1,   -0.1,   -0.1,   -0.1,   -0.1,    0]'; % zl
        
        maxSec = [1, 0.9, 0.7,  0.5,    0.3,    0,      0;  % xu
            1,      0.9,    0.7,    0.5,    0.3,    0,      0;  % xl
            0.03,   0.1,    0.1,    0.1,    0.1,    0.1,    0;  % zu
            0,      0,      0,      0,      0,      0,      0]'; % zl
        
        % Use to hold sections at base/standard throughout simulation
        %     [minSec,maxSec] = deal(nan(BezierControlPoints + 1,4));
        
        [minControlPoints,~] = size(minSec(2:end,:));
        [maxControlPoints,~] = size(maxSec(2:end,:));
        
        if any([minControlPoints,maxControlPoints] ~=  BezierControlPoints)
            
            error("Number of Bezier control points not consistent with options definition")
        end
        
        minSec = repmat(minSec, 1, 1, n + 1);
        maxSec = repmat(maxSec, 1, 1, n + 1);
        
        options.nBez = BezierControlPoints + 1;
        
        sectionStr = "Bezier";
        sectionCons = "~";
        sectionTrans = "~";
        
        foilDefs = {...
            "Variables",    "VarMin",   "VarMax",   "Conditions",	"Transformations",	"Optimise/Hold"
            sectionStr,    minSec,     maxSec,     sectionCons,	sectionTrans,       "Optimise"};
        
    case "BezierTC"
        
        BezierControlPoints = options.BezierControlPoints;
        
        % Control point min/max coordinates
%         minSec = [1, 0.1,   0.1,    0.1,    0.1,    0,      0;  % xc
%             1,      0.1,    0.1,    0.1,    0.1,    0,      0;  % xt
%             -0.03,  -0.1,   -0.1,   -0.05,  -0.05,   0,      0;  % zc
%             0,     -0.05,  -0.05,  -0.05,  -0.05,   0.05,   0]'; % zt
%         
%         maxSec = [1, 0.9,   0.9,    0.9,    0.9,    0,      0;  % xc
%             1,      0.9,    0.9,    0.9,    0.9,    0,      0;  % xt
%             0.03,   0.1,    0.2,    0.2,    0.2,    0.05,   0;  % zc
%             0.03,   0.1,    0.2,    0.2,    0.2,    0.1,    0]'; % zt
        
%         minSec = [1, 0.7,   0.5,    0.3,    0.1,    0,      0;  % xc
%             1,      0.7,    0.5,    0.3,    0.1,    0,      0;  % xt
%             -0.03,  -0.1,   -0.1,   -0.05,  -0.05,   0,      0;  % zc
%             0,     -0.05,  -0.05,  -0.05,  -0.05,   0.05,   0]'; % zt
%         
%         maxSec = [1, 0.9,   0.7,    0.5,    0.3,    0,      0;  % xc
%             1,      0.9,    0.7,    0.5,    0.3,    0,      0;  % xt
%             0.03,   0.1,    0.2,    0.2,    0.2,    0.05,   0;  % zc
%             0.03,   0.1,    0.2,    0.2,    0.2,    0.1,    0]'; % zt
        
        %% NEW TEST
        minSec = [1, 0.7,   0.5,    0.3,    0.1,    0,      0;  % xc
            1,      0.7,    0.5,    0.3,    0.1,    0,      0;  % xt
            -0.03, -0.1,   -0.1,   -0.05,  -0.05,   0,      0;  % zc
            0,     -0.05,  -0.05,  -0.05,  -0.05,   0.01,   0]'; % zt
        
        maxSec = [1, 0.9,   0.7,    0.5,    0.3,    0,      0;  % xc
            1,      0.9,    0.7,    0.5,    0.3,    0,      0;  % xt
            0.03,   0.1,    0.2,    0.2,    0.2,    0.01,   0;  % zc
            0.03,   0.1,    0.2,    0.2,    0.1,    0.03,   0]'; % zt
        
        %%
        % Use to hold sections at base/standard throughout simulation
        %     [minSec,maxSec] = deal(nan(BezierControlPoints + 1,4));
        
        [minControlPoints,~] = size(minSec(2:end,:));
        [maxControlPoints,~] = size(maxSec(2:end,:));
        
        if any([minControlPoints,maxControlPoints] ~=  BezierControlPoints)
            
            error("Number of Bezier control points not consistent with options definition")
        end
        
        minSec = repmat(minSec, 1, 1, n + 1);
        maxSec = repmat(maxSec, 1, 1, n + 1);
        
        options.nBez = BezierControlPoints + 1;
        
        sectionStr = "Bezier";
        sectionCons = "~";
        sectionTrans = "~";
        
        foilDefs = {...
            "Variables",    "VarMin",   "VarMax",   "Conditions",	"Transformations",	"Optimise/Hold"
            sectionStr,    minSec,     maxSec,     sectionCons,	sectionTrans,       "Optimise"};
        
    case "Preloaded"
        
        % Define sections be pre-loaded data files
        minSec = 1;
        
        % Load coordinates of 2D aerofoil sections into matrices within cell,
        % max defined by number of stored data files
        [options.foilData,maxSec] = getaerofoilsecdata();
        
        sectionStr = "Section";
        sectionCons = "Floor";
        sectionTrans = "~";
        
        foilDefs = {...
            "Variables",    "VarMin",   "VarMax",   "Conditions",	"Transformations",	"Optimise/Hold"
            sectionStr,    minSec,     maxSec,     sectionCons,	sectionTrans,       "Optimise"};
        
    otherwise
        error('No method available for user-specified aerofoil creation method %s', aerofoilMethod)
end

foilDefs = optorhold(foilDefs);