function [Cp,Mach,P] = tangentempirical(del,flow,conical)

alpha=flow.alpha; Minf=flow.Minf; gamma=flow.gamma; Pinf = flow.Pinf;

if alpha <= 15 * pi/180 && ~conical
    
    Kw = (gamma+1)/2;
    a = Kw*Minf*sin(del);
    Mns = a+exp(-a)/2;
else
    Kc = (2*(gamma+1))/(gamma+3);
    a = Kc*Minf*sin(del);
    Mns = a+exp(-a);
end

% Cp = 2*(sin(del).^2).*(1-(((gamma-1).*(Mns.^2)+2)./(4.*(gamma+1).*(Mns.^2)))).^-1;
Cp = 2*(sin(del).^2).*(1-(((gamma-1).*(Mns.^2)+2)./(4.*(gamma+1).*(Mns.^2))));

P = (2*Cp .* gamma .* Pinf .* Minf^2) + Pinf;

Mach = 

con = isnan(Cp);

Cp(con) = 0;
Mach(con) = 0;