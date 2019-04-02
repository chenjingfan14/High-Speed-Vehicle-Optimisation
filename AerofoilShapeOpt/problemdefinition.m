%% Initialise aerofoil shape optimisation
clear all
close all
clc

restoredefaultpath

addpath(genpath('AerofoilShapeOpt'))
addpath(genpath('Optimisation'))
addpath(genpath('VehicleGen'))
addpath(genpath('Functions'))

resultPath = [pwd '\Results'];

baseline = false;
cluster = false;

%%
% Number of processors to be used (up to 4 on desktop)
nProc = 1;
costFun = @aerofoilcompare; % Cost function caller

%% Simulation options

% Number of cost functions
nFun = 1;

% Aerofoil creation method: Preloaded, Bezier, BP3434
aerofoilMethod = 'BezierTC';

% Number of Bezier control points to use if that is the chosen method
BezierControlPoints = 6;

% Number of wing partitions (for this optimisation it is 0 ie 2D aerofoil)
wingPartitions = 0;

tolerance = 8e-4;

%%

options.CostFunctions = nFun;

% If number of processors has been entered and if that value is > 1, create
% parallel loop
if exist('nProc','var') && nProc > 1
    options.Parallel = true;
    
    pool = gcp('nocreate');

    if isempty(pool)
        parpool('local',nProc);
    end
    
else % Running on one processor
    options.Parallel = false;
end

options.AerofoilMethod = aerofoilMethod;
options.WingPartitions = wingPartitions;
options.HardTransform = true;
options.Tolerance = tolerance;
options.Baseline = baseline;
options.Cluster = cluster;

% Failsafe incase forget to change above when running on cluster
if exist('nProc','var') && nProc > 4
    options.Cluster = true;
end

if contains(aerofoilMethod,"Bezier")
    
    options.BezierControlPoints = BezierControlPoints;
end

%% Initialise optimisation variables

%% Part Definitions
variCons = {...
    "Variables", "VarMin", "VarMax", "Conditions", "Transformations", "Optimise/Hold"};

foilDefs = define_aerofoil(options);
variCons = [variCons; foilDefs];

%% Baseline
% If there is a baseline configuration to be bettered, create here and
% provide baseline parameters/characteristics
if baseline
    
    options.Base = baselinefun(variCons,options);
end

%% Output
% variCons = standardvariables(variCons,options);

% Streamlines above cell to indices instead of names
% Total number of variables
[cond,varArray,varMin,varMax,nVar] = translateOpt(variCons);