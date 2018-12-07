clear all
close all
clc

addpath(genpath('MultiSwarm'))
addpath(genpath('VehicleGen'))
addpath(genpath('AeroPrediction'))

%%
% Number of processors to be used (up to 4 on desktop)
nProc = 1;
costFun = @aeroprediction; % Cost function caller
options = simOptions(nProc);

Wing = options.Wing;
Aft = options.Aft;
Fore = options.Fore;
Nose = options.Nose;
Control = options.Control;

parallel = options.Parallel;
Bezier = options.Bezier;
baseline = options.Baseline;
cluster = options.Cluster;

% Must be called with options, additional inputs can be (in this order)
% arrays of: angle of attack, Mach, altitude, control surface deflection,
% although pre-determined values for these can be set in flowparameters.
% For example: use pre-determined AoA, 2 Mach numbers, one altitude, no
% control - flowparameters(options,[],[3,4],[10000])
flow = flowparameters(options);

% Load lookup tables for shock-expansion and Prandtl Meyer expansion
load('thetaBetaCurves.mat');
Mrange = [1:0.0001:10,10.1:0.1:100];
PrandtlMeyer = prandtlmeyerlookup(Mrange,flow);

% Cell containing variable names and their subsequent conditions to be
% applied, so far these conditions are as follows:
% minimum value (mini) - if value is below those thresholds, set to zero
% floor - round value down to nearest integer
% if all - if any values = 0, set all to zero

% Number of wing partitions
n = 3;

%% Part Definitions
variCons = {...
    "Variables",            "VarMin",   "VarMax",   "Num Of",       "Conditions"                "Transformations"};

wingDefs = {...
    "Dihedral",             0,          20,         "~",            "~",                        "~";
    "Chord",                0.1,        1,          n+1,            "< Previous",               ".*AftLength";
    "LESweep",              0,          80          n,              "~"                         "~";
    "Semispan",             0,          5,          n,              ["Minimum 0.5" "sum > 2"],  "~";
    "SectionDefinition",    [],         [],         "~",            "~",                        "~";
    "xOffset",             -0.25,       0.5,        "~",            "~",                        ".*AftLength";
    "zOffset",             -0.5,        0,          "~",            "~"                         ".*AftHalfHeight"};

aftDefs = {...
    "UpperLength",          "~",        "~",        "~",            "~",                        "~";
    "yUpperRad",            "~",        "~",        "~",            "~",                        "~";
    "yBotRatio",            "~",        "~",        "~",            "~",                        "~";
    "zUpperRad",            "~",        "~",        "~",            "~",                        "~";
    "SideLength",           "~",        "~",        "~",            "~",                        "~";
    "zLowerRad",            "~",        "~",        "~",            "~",                        "~";
    "AftLength",            "~",        "~",        "~",            "~",                        "~"};

foreDefs = {...
    "ForeLength",           "~",        "~",        "~",            "~",                        "~"};

noseDefs = {...
    "NoseRad",              "~",        "~",        "~",            "~",                        "~";
    "NoseLength",           "~",        "~",        "~",            "~",                        ".*NoseRad";
    "zNoseOffset",          "~",        "~",        "~",            "~",                        ".*AftHalfHeight"};

controlDefs = {...
    "ControlSpan",          "~",        "~",        2,              "> Previous",               "~";
    "ControlChord",         "~",        "~",        "~",            "~",                        "~"};

%% 

if Wing
    
    % 2D Aerofoil Section Defintion
    
    % Define sections as Bezier curves
    if Bezier
        
        % Control point min/max coordinates
        minSec = [1, 0.7, 0.5,  0.3,    0.1,    0,      0;  % xu
            1,  0.7,    0.5,    0.3,    0.1,    0,      0;  % xl
            0,  0.015,  0.02,   0.05,   0.02,   0.05    0;  % zu
            0, -0.035, -0.04,  -0.07,  -0.04,  -0.05,   0]; % zl
        
        maxSec = [1, 0.9, 0.7,  0.5,    0.3,    0,      0;  % xu
            1,  0.9,    0.7,    0.5,    0.3,    0,      0;  % xl
            0,  0.035,  0.04,   0.07,   0.04,   0.05,   0;  % zu
            0, -0.015, -0.02,  -0.05,  -0.02,  -0.05,   0]; % zl
        
        % Transform matrix to single array
        minSec = reshape(minSec',1,[]);
        maxSec = reshape(maxSec',1,[]);
        
        foilData = length(minSec);
        
        SectionDefinition = {"Bezier", minSec, maxSec, n+1, "~", "~"};
        
    else % Define sections be pre-loaded data files
        
        minSec = 1;
        
        % Load coordinates of 2D aerofoil sections into matrices within cell,
        % max defined by number of stored data files
        [foilData,maxSec] = getaerofoilsecdata();
        
        SectionDefinition = {"Section", minSec, maxSec, n+1, "Floor", "~"};
        
    end
    
    % Find which set of variables correspond to 2D aerofoil sections
    for i=size(wingDefs,1):-1:1
        secDef(i) = wingDefs{i,1} == "SectionDefinition";
    end

    % And Insert definition
    wingDefs(secDef,:) = SectionDefinition;
    
    variCons = [variCons; wingDefs];
    
end

if Aft
    variCons = [variCons; aftDefs];
end

if Fore
    variCons = [variCons; foreDefs];
end

if Nose
    variCons = [variCons; noseDefs];
end

if Control
    variCons = [variCons; controlDefs]; 
end

% Streamlines above cell to indices instead of names
% Total number of variables
[cond,varArray,varMin,varMax,nVar] = translateOpt(variCons);

%%

% If there is a baseline configuration to be bettered, create here and
% provide baseline parameters/characteristics
if baseline
    options.base = baselinefun(flow,options,thetaBetaM,maxThetaBetaM,PrandtlMeyer);
end

standard = isnan(varMin)';
[varMin,varMax] = standardvariables(standard,n,options,varMin,varMax,nVar,varArray);