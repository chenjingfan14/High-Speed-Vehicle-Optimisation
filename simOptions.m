function options = simOptions(nProc)
%% Include what in configurations
Wing = true; Aft = true; Fore = true; Nose = true;

options.Wing = Wing;
options.Aft = Aft;
options.Fore = Fore;
options.Nose = Nose;

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

% Use Bezier splines for 2D aerofoil definition, else use preloaded data
% files
options.Bezier = true;

% Use hard coded transforms (required for cluster), else uses versatile
% (sym engine required) transform function
options.Hardtransform = true;

% Leave false, needs validated
options.Shielding = false;

% Include viscous effects in aerodynamic prediction
options.Viscous = true;

% Include control surfaces as design variables
options.Control = false;

% Initial baseline configuration to be analysed and used as reference
options.Baseline = true;

% Running simulation on computing cluster (ie. Buckethead)?
options.Cluster = false;

% Failsafe incase forget to change above when running on cluster
if exist('nProc','var') && nProc > 4
    options.Cluster = true;
end