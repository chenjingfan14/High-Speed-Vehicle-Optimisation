function [Cp] = vandykeunified(del,flow)

Minf = flow.Minf; gamma = flow.gamma;

H = del * (Minf^2 - 1)^0.5;

% Compression
Cp = (del.^2)*(((gamma+1)/2) + (((gamma+1)/2)^2 + (4/(H^2)))^0.5);

% Expansion (No LE shock)
Cp = (del.^2) * (2/(gamma*H^2)) * ((1 - ((gamma-1)/2) * H)^((2*gamma)/(gamma-1)) - 1);