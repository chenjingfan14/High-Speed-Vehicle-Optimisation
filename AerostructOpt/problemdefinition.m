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

%% Simulation options

% "PSO" or "GA"
method = "PSO";

% Number of processors to be used (up to 4 on desktop)
nProc = 1;
costFun = @particlecreator; % Cost function caller

% Number of cost functions
nFun = 2;

% If opt cost function value = max(f(x)) (rather than min) then can use
% this to invert CF values for display purposes
options.Inv = false(1,nFun);
% options.Neg = false(1,nFun);

options.Neg = logical([1 0]);

options.CostFunctions = nFun;
options.CostFunctionLabels = {'C_L','C_D'};

% Tolerance required to end optimisation (single-objective only)
options.Tolerance = 0;

%% Optimisation inclusions and methods
options.Wing = true;
options.Aft = true;
options.Fore = true;
options.Nose = true;
options.Control = false;
options.Quad = false;

% Prediction method mixer:
    % Impact:   1 - Modified Newtonian
    %           2 - Modified Newtonian + Prandtl-Meyer
    %           3 - Oblique Shock + Prandtl-Meyer
    %           4 - Tangent Wedge/Cone
    % Shadow:   1 - Newtonian/Base Pressure
    %           2 - Prandtl-Meyer
    %
    % Bodypart: whether or not part should be treated as body (true) or
    % lifting surface
aeroMethod = {["aerofoil","wing","tail","aerofoil l","aerofoil u"], 3,  2,  false;
               "aftbody",                                           3,  2,  true;
               "forebody",                                          3,  2,  true;
               "nose",                                              3,  2,  true;
               "test",                                              3,  2,  true;
               "default",                                           3,  2,  true};

options.AeroMethod = aeroMethod;

% Aerofoil creation method: PARSEC, Bezier, BP3434, Preloaded
aerofoilMethod = "BezierTC";

% Number of Bezier control points to use if that is the chosen method
BezierControlPoints = 6;

% Number of wing partitions
n = 2;

% Chordwise wing discretisation
xPanels = 20;
discMethod = "Cosine";
          
% Use hard coded transforms (required for cluster), or versatile
% (sym engine required) transform function
options.HardTransform = true;

% Include structure (not yet complete)
options.Structure = false;

% Leave false, needs validated
options.Shielding = false;

% Include viscous effects in aerodynamic prediction
options.Viscous = false;

% "Euclids Hat" "Wendlands" "Thin Plate Spline". Default is Volume Spline
rbfmethod = "Wendlands";

% Initial baseline configuration to be analysed and used as reference
options.Baseline = true;

%% Flow parameters

flow.Alpha = (-4:4:24)';
flow.Alpha = [22.8 23.3]';
flow.Delta = 0;
flow.Control = options.Control;

% Nose-body Mach 4.63
% flow.Mach = 4.63;
% flow.Altitude = 32850;

% X-34 Mach 6
% flow.Mach = 6;
% flow.Altitude = 27975;

% X-34 Mach 3 (2?)
flow.Mach = 2;
flow.Altitude = 20000;
flow.Altitude = 13530;

flow = flowparameters(flow);


%% Creating simulation parameters based on chosen options

% If number of processors has been entered and if that value is > 1, create
% parallel loop
if exist('nProc','var') && nProc > 1
    
    options.Parallel = true;
    
    % Turn on when running on cluster to avoid graphical output
    if nProc > 4
    
        options.Cluster = true;
    end
        
    pool = gcp('nocreate');

    if isempty(pool)
        parpool('local',nProc);
    end
    
else % Running on one processor
    options.Parallel = false;
end

if options.Structure
    
    opts.SYM = true;
    opts.POSDEF = true;
    options.StructOpts = opts;
end

if options.Viscous
    
    options.RBFfun = rbffunctions(rbfmethod);
end

%% Initialise optimisation variables
aft = options.Aft;
fore = options.Fore;
nose = options.Nose;

% Load lookup tables for shock-expansion and Prandtl Meyer expansion
% Mrange = [1:0.0001:10,10.1:0.1:100];
% options.PrandtlMeyer = prandtlmeyerlookup(Mrange,flow);
% options.PrandtlMeyer = [];

% Can insert own (Mrange,betaRange) vals into function
[~, options.MaxThetaBetaM] = thetabetamachcurves();

options.pmFun = @(M1,gamma) ((gamma + 1)/(gamma - 1)).^0.5 * atan((((gamma - 1)/(gamma + 1)).*(M1.^2 - 1)).^0.5) - atan(((M1.^2) - 1).^0.5);

%% Part Definitions
variCons = {...
    "Variables", "VarMin", "VarMax", "Conditions", "Transformations", "Optimise/Hold"};

if options.Wing
    
    X = (0:xPanels)';
    
    switch discMethod

        case "Linear"

            chordDisc = X/max(X);

        case "Cosine"

            chordDisc = 0.5*(1-cos((X*pi)/max(X)));

        case "HalfCosine"

            chordDisc = 1-cos((X*(pi/2))/max(X));
    end
    
    if contains(aerofoilMethod,'Bezier')
        
        options.BezierControlPoints = BezierControlPoints;
    end
    
    options.WingPartitions = n;
    options.AerofoilMethod = aerofoilMethod;
    options.ChordDisc = chordDisc;
    
    wingDefs = define_wing(options);
    [foilDefs,options] = define_aerofoil(options);
    variCons = [variCons; wingDefs; foilDefs];
else
    
    foilData = NaN;
end

if any([aft,fore,nose])
    
    bodyDefs = define_body(aft,fore,nose);
    variCons = [variCons; bodyDefs];    
end

if options.Control
    
    controlDefs = define_control();
    variCons = [variCons; controlDefs];
end

options.Flow = flow;

%% Baseline

% If there is a baseline configuration to be bettered, create here and
% provide baseline parameters/characteristics
if options.Baseline
    
    options.Base = baselinefun(variCons,options);
end

%% Output
variCons = standardvariables(variCons,options);

% Streamlines above cell to indices instead of names
% Total number of variables
[cond,varArray,varMin,varMax,nVar] = translateOpt(variCons);