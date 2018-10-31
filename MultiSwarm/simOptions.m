function options = simOptions()

% For parallel processing set as true, otherwise false
options.parallel = false;
options.shielding = false;
options.Bezier = false;

% Viscous capability still being created, leave false
options.viscous = false;