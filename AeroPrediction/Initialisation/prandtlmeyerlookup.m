function PrandtlMeyer = prandtlmeyerlookup(Mrange,flow)

gamma = flow.gamma;

frac = (gamma+1)/(gamma-1);

vM = (frac^0.5)*atan(((1/frac)*((Mrange.^2)-1)).^0.5) - atan(((Mrange.^2)-1).^0.5);

PrandtlMeyer = [Mrange;vM]';

% plot(Mrange,vM)