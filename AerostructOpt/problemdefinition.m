%% Initialise aero-structural optimisation problem, paths etc

clear all
close all
clc

restoredefaultpath

addpath(genpath('AerostructOpt'))
addpath(genpath('Optimisation'))
addpath(genpath('VehicleGen'))
addpath(genpath('Aerodynamics'))
addpath(genpath('Structures'))
addpath(genpath('Functions'))

resultPath = [pwd '\Results'];

%%
% Number of processors to be used (up to 4 on desktop)
nProc = 1;
costFun = @particlecreator; % Cost function caller

%% Simulation options

% Number of cost functions
nFun = 1;

%% Include what in optimisation
options.Wing = true;
options.Aft = true;
options.Fore = true;
options.Nose = true;
options.Control = false;

% Aerofoil creation method: Preloaded, Bezier, BP3434
aerofoilMethod = "BezierTC";

% Number of Bezier control points to use if that is the chosen method
BezierControlPoints = 6;

% Number of wing partitions
n = 3;

% Chordwise wing discretisation
xPanels = 50;
discMethod = "Cosine";

X = (0:xPanels)';

switch discMethod
    
    case "Linear"
        
        chordDisc = X/max(X);
        
    case "Cosine"
        
        chordDisc = 0.5*(1-cos((X*pi)/max(X)));
        
    case "HalfCosine"
        
        chordDisc = 1-cos((X*(pi/2))/max(X));
end

%% Flow parameters

flow.Alpha = (-4:4:24)';
flow.Mach = 4.63;
flow.Altitude = 32850;
flow.Delta = 0;
flow.Control = options.Control;

flow = flowparameters(flow);

%%

options.Tolerance = 0;

% Use hard coded transforms (required for cluster), else uses versatile
% (sym engine required) transform function
options.HardTransform = true;

% Include structure (not yet complete)
options.Structure = false;

% Leave false, needs validated
options.Shielding = false;

% Include viscous effects in aerodynamic prediction
options.Viscous = false;

% % Include control surfaces as design variables
% options.Control = control;

% Initial baseline configuration to be analysed and used as reference
options.Baseline = true;

% Running simulation on computing cluster (ie. Buckethead)?
options.Cluster = false;

%%

options.CostFunctions = nFun;

% If number of processors has been entered and if that value is > 1, create
% parallel loop
if exist('nProc','var') && nProc > 1
    
    parallel = true;
    
    pool = gcp('nocreate');

    if isempty(pool)
        parpool('local',nProc);
    end
    
else % Running on one processor
    parallel = false;
end

options.Parallel = parallel;

if options.Wing
    
    options.WingPartitions = n;
    options.AerofoilMethod = aerofoilMethod;
end

if contains(aerofoilMethod,'Bezier')
    
    options.BezierControlPoints = BezierControlPoints;
end

% Failsafe incase forget to change above when running on cluster
if exist('nProc','var') && nProc > 4
    
    options.Cluster = true;
end

%% Initialise optimisation variables
wing = options.Wing;
aft = options.Aft;
fore = options.Fore;
nose = options.Nose;
control = options.Control;

baseline = options.Baseline;
cluster = options.Cluster;

% Load lookup tables for shock-expansion and Prandtl Meyer expansion
load('thetaBetaCurves.mat');
Mrange = [1:0.0001:10,10.1:0.1:100];
options.PrandtlMeyer = prandtlmeyerlookup(Mrange,flow);
options.MaxThetaBetaM = maxThetaBetaM;
options.ThetaBetaM = thetaBetaM;

% Cell containing variable names and their subsequent conditions to be
% applied, so far these conditions are as follows:
% minimum value (mini) - if value is below those thresholds, set to zero
% floor - round value down to nearest integer
% if all - if any values = 0, set all to zero

%% Part Definitions
variCons = {...
    "Variables", "VarMin", "VarMax", "Conditions", "Transformations", "Optimise/Hold"};

if wing
    
    wingDefs = define_wing(options);
    [foilDefs,opotions] = define_aerofoil(options);
    variCons = [variCons; wingDefs; foilDefs];
else
    
    foilData = NaN;
end

if any([aft,fore,nose])
    
    bodyDefs = define_body(aft,fore,nose);
    variCons = [variCons; bodyDefs];    
end

if control
    
    controlDefs = define_control();
    variCons = [variCons; controlDefs];
end

options.Flow = flow;
options.ChordDisc = chordDisc;

%% Baseline

% If there is a baseline configuration to be bettered, create here and
% provide baseline parameters/characteristics
if baseline
    
    options.Base = baselinefun(variCons,options);
end

%% Output
variCons = standardvariables(variCons,options);

% Streamlines above cell to indices instead of names
% Total number of variables
[cond,varArray,varMin,varMax,nVar] = translateOpt(variCons);