%% Initialise problem, paths, 

clear all
close all
clc

addpath(genpath('MultiSwarm'))
addpath(genpath('VehicleGen'))
addpath(genpath('AeroPrediction'))
addpath(genpath('Structures'))

%%
% Number of processors to be used (up to 4 on desktop)
nProc = 1;
costFun = @aeroprediction; % Cost function caller
options = simOptions(nProc);
% Number of decision variables (cost function values)
nFun = options.CostFunctions;

wing = options.Wing;
aft = options.Aft;
fore = options.Fore;
nose = options.Nose;
control = options.Control;
n = options.WingPartitions;

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

%% Part Definitions
variCons = {...
    "Variables", "VarMin", "VarMax", "Conditions", "Transformations", "Optimise/Hold"};

%% 

if wing
    
    [wingDefs,foilData] = define_wing(options);
    variCons = [variCons; wingDefs];
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

%%

% If there is a baseline configuration to be bettered, create here and
% provide baseline parameters/characteristics
if baseline
    
    options.Base = baselinefun(variCons,flow,options,thetaBetaM,maxThetaBetaM,PrandtlMeyer);
end

variCons = standardvariables(variCons,options);

% Streamlines above cell to indices instead of names
% Total number of variables
[cond,varArray,varMin,varMax,nVar] = translateOpt(variCons);