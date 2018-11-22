function options = simOptions(nProc)

if nProc > 1
    options.parallel = true;
else
    % For parallel processing set as true, otherwise false
    options.parallel = false;
end

% Use Bezier splines for 2D aerofoil definition, else use preloaded data
% files
options.Bezier = true;

% Use hard coded conditions (required for cluster), else uses versatile
% conditioning function
options.hardconditioning = true;

% Leave false, needs validated
options.shielding = false;

% Viscous capability still being created, leave false
options.viscous = true;

% Include control surfaces in optimisation
options.control = false;

% Initial baseline configuration to be analysed and used as reference vals
options.baseline = true;