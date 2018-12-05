function options = simoptions(nProc)

% If number of processors has been entered and if that value is > 1, create
% parallel loop
if exist('nProc','var') && nProc > 1
    options.parallel = true;
    
    pool = gcp('nocreate');

    if isempty(pool)
        parpool('local',nProc);
    end
    
else % Running on one processor
    options.parallel = false;
end

% Use Bezier splines for 2D aerofoil definition, else use preloaded data
% files
options.Bezier = true;

% Use hard coded transforms (required for cluster), else uses versatile
% (sym engine required) transform function
options.hardtransform = true;

% Leave false, needs validated
options.shielding = false;

% Include viscous effects in aerodynamic prediction
options.viscous = true;

% Include control surfaces as design variables
options.control = false;

% Initial baseline configuration to be analysed and used as reference
options.baseline = true;

% Running simulation on computing cluster (ie. Buckethead)?
options.cluster = false;

% Failsafe incase forget to change above when running on cluster
if exist('nProc','var') && nProc > 4
    options.cluster = true;
end