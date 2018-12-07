clear all
close all
clc

addpath(genpath('MultiSwarm'))
addpath(genpath('VehicleGen'))
addpath(genpath('AeroPrediction'))

% Number of processors to be used (up to 4 on desktop)
nProc = 1;
costFun = @aeroprediction; % Cost function caller
options = simOptions(nProc);

parallel = options.parallel;
Bezier = options.Bezier;
control = options.control;
baseline = options.baseline;
cluster = options.cluster;

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

variCons = {"Variables",    "Num Of",   "Conditions"    'Transformations';...
    "Dihedral",             "~",            "~",                        '~';...
    "Chord",                n+1,            "< Previous",               '.*AftLength';...
    "LESweep",              n,              "~"                         '~';...
    "Semispan",             n,              ["Minimum 0.5" "sum > 2"],  '~';...            
    "SectionDefinition",    "~",            "~",                        '~';...
    "xOffset",              "~",            "~",                        '.*AftLength';...
    "zOffset",              "~",            "~"                         '.*AftHalfHeight';...
    "UpperLength",          "~",            "~",                        '~';...
    "yUpperRad",            "~",            "~",                        '~';...
    "yBotRatio",            "~",            "~",                        '~';...
    "zUpperRad",            "~",            "~",                        '~';...
    "SideLength",           "~",            "~",                        '~';...
    "zLowerRad",            "~",            "~",                        '~';...
    "AftLength",            "~",            "~",                        '~';...
    "NoseRad",              "~",            "~",                        '~';...
    "NoseLength",           "~",            "~",                        '.*NoseRad';...
    "zNoseOffset",          "~",            "~",                        '.*AftHalfHeight';...
    "ForeLength",           "~",            "~",                        '~'};

%% Minimum/Maximum Variable Definitions
% Corresponds to above. Current setup has all partitions with equal min/max
% values for partitions dependent variables
minChord = 0.1;
minLESweep = 0;
minSemispan = 0;

minChord = zeros(1,n+1) + minChord;
minLESweep = zeros(1,n) + minLESweep;
minSemispan = zeros(1,n) + minSemispan;

maxChord = 1;
maxLESweep = 80;
maxSemispan = 5;

maxChord = zeros(1,n+1) + maxChord;
maxLESweep = zeros(1,n) + maxLESweep;
maxSemispan = zeros(1,n) + maxSemispan;

%% 2D Aerofoil Section Defintion

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

    SectionDefinition = {"Bezier", foilData*(n+1), "~", '~'};

else % Define sections be pre-loaded data files
    
    minSec = 1;
    
    % Load coordinates of 2D aerofoil sections into matrices within cell,
    % max defined by number of stored data files
    [foilData,maxSec] = getaerofoilsecdata();
    
    SectionDefinition = {"Section", n+1, "Floor", '~'};

end

% Find which set of variables correspond to 2D aerofoil sections
for i=size(variCons,1):-1:1
    secDef(i) = variCons{i,1} == "SectionDefinition";
end

% And Insert definition
variCons(secDef,:) = SectionDefinition;

% Replicate min/max section values for number of sections ie. number of
% partitions+1
minSec = repmat(minSec,1,n+1);
maxSec = repmat(maxSec,1,n+1);

% Minimum and maximum variable values. Corresponds to variCons
%% Full Configuration (Wing & Body)
% varMin = [0, 0.5,0.1,0.1,0.1, 0,0,0, 2,0,0, 1,1,1,1, 0,-0.5,... % Wing
%     0,0.1,0.1, 0.1,0,0.1, 4, 0,0,-0.5, 1]; % Body
% varMax = [20, 1,1,1,1, 45,45,45, 5,2,2, 33,33,33,33, 0.5,0,...
%     1,1,1, 1,1,1, 10, 0.5,1,0, 5];

%% Wing Only
varMin = [0, minChord, minLESweep, minSemispan, minSec, -0.25,-0.5,... % Wing
    NaN,NaN,NaN, NaN,NaN,NaN, NaN, NaN,NaN,NaN, NaN]; % Body
varMax = [20, maxChord, maxLESweep, maxSemispan, maxSec, 0.5,0,...
    NaN,NaN,NaN, NaN,NaN,NaN, NaN, NaN,NaN,NaN, NaN];

%% Control Parameters
if control
    
    minControl = [0.1,0.2,0.5];
    maxControl = [0.9,0.99,0.9];
    
    variCons(end+1,:) = {"ControlSpan", 2, "> Previous", '~'};
    variCons(end+1,:) = {"ControlChord", "~", "~", '~'};
    
    varMin = [varMin, minControl];
    varMax = [varMax, maxControl];
    
end

% Streamlines above cell to indices instead of names
% Total number of variables
[cond,varArray,nVar] = translate(variCons);

%%

% If there is a baseline configuration to be bettered, create here and
% provide baseline parameters/characteristics
if baseline
    options.base = baselinefun(flow,options,thetaBetaM,maxThetaBetaM,PrandtlMeyer);
end

if exist('varMin','var')
    standard = isnan(varMin);
    [varMin,varMax] = standardvariables(standard,n,options,varMin,varMax);
else
    standard = true(1,nVar);
    [varMin,varMax] = standardvariables(standard,n,options);
end