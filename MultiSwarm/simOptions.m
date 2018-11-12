function options = simOptions()

% For parallel processing set as true, otherwise false
options.parallel = false;

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