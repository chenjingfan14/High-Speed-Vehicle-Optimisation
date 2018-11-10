function Cdf = friction(points,flow)

%% Eckert's Reference Temperature method
% Iteration used in Thermoelastic Formulation of a Hypersonic Vehicle
% Control Surface for Control-Oriented Simulation

% Constant 0.7 in above ref
Pr = flow.Pr;
Me = flow.Mach;
Te = flow.Tinf;
gamma = flow.gamma;

r = Pr.^(1/3);

Tt = Te*(1 + (gamma - 1)*((Me^2)/2));

Tr = r*(Tt - Te) + Te;

TRef = Te + 0.5*(Tw - Te) + 0.22*(Tr - Te);

% Ideal gas law
rhoRef = 0;
% Sutherland's law
muRef = 0;

ReRefx = (rhoRef*Ve*x)/muRef;

% Accurate up to 10^9 according to above ref
cf = 0.37/(log(ReRefx).^2.584);