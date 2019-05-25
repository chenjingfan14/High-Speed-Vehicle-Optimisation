function Cp = basepressure(Minf,gamma)

% Cp = 0;
% Cp = -1./(Minf.^2);
Cp = (2./(gamma * Minf.^2)) .* ((2/(gamma + 1))^1.4 .* (1./Minf).^2.8 .* ((2 * gamma * (Minf.^2) - (gamma - 1))/(gamma + 1)) - 1);
