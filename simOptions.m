function options = simOptions(nProc)
%% Simulation options

% Number of cost functions
nFun = 1;

%% Include what in optimisation
options.Wing = true;
options.Aft = true;
options.Fore = true;
options.Nose = true;
options.Control = false;

% Number of wing partitions
wingPartitions = 3;

% Use Bezier splines for 2D aerofoil definition, else use preloaded data
% files
options.Bezier = true;
% Number of Bezier control points
BezierControlPoints = 6;

%% Flow parameters

options.Flow.Alpha = (-4:4:24)';
options.Flow.Mach = 4.63;
options.Flow.Altitude = 32850;
options.Flow.Delta = 0;

%%

% Use hard coded transforms (required for cluster), else uses versatile
% (sym engine required) transform function
options.Hardtransform = true;

% Include structure (not yet complete)
options.Structure = false;

% Leave false, needs validated
options.Shielding = false;

% Include viscous effects in aerodynamic prediction
options.Viscous = true;

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
    options.Parallel = true;
    
    pool = gcp('nocreate');

    if isempty(pool)
        parpool('local',nProc);
    end
    
else % Running on one processor
    options.Parallel = false;
end

if options.Wing
    options.WingPartitions = wingPartitions;
end

if options.Bezier
    options.BezierControlPoints = BezierControlPoints;
end

% Failsafe incase forget to change above when running on cluster
if exist('nProc','var') && nProc > 4
    options.Cluster = true;
end