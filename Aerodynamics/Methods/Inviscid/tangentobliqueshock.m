function [Cp,Mach,P,method] = tangentobliqueshock(del,ID,method,flow,maxThetaBetaM,conical)

M1 = flow.Minf;
gamma = flow.gamma;
Pinf = flow.Pinf;

[Cp,Mach,P] = deal(zeros(size(del)));

target = tan(del);

ub = halfspace(M1,maxThetaBetaM(:,1),maxThetaBetaM(:,3));

tbm = @(B,M1,gamma) 2*cot(B) .* ((M1.^2) .* (sin(B).^2) - 1)./((M1.^2) .* (gamma + cos(2*B)) + 2);
tbmbi = @(B) tbm(B,M1,gamma);

[beta, ~, flag] = bisection(tbmbi,0,ub,target(:));

% All panel inclinations to flow for attached shock must be less than
% maximum allowable otherwise switch infringing panels to newtonian flow
con = ~isnan(beta);

beta = beta(con);

thetaO = del(con);
thetaN = del(~con);

M = zeros(sum(con),1);

if conical
    % Development of an Aerodynamics Code for the Optimisation of 
    % Hypersonic Vehicles - gives M2 > M1
    % Mix between On hypersonic flow past unyawed cone & Approximate
    % solutions for supersonic flow over wedges and cones

    tau = asin(sin(thetaO) .* (((gamma+1)/2) + (1./((M1*sin(thetaO)).^2))).^0.5);

    % numer = (Minf^2) .* (cos(tau).^2) .* (1+2*(tau-theta));
    % denom = 1 + ((gamma-1)/2) * (Minf^2) * ((sin(tau).^2) - 2*((tau-theta).^2) .* (cos(theta).^2));
    % M = (numer./denom).^0.5;
    
    % Exact Taylor-Maccoll solution?
    for i = 1:length(tau)
        
        [~,M(i),~] = solvecone(tau(i),M1,gamma);
    end
    
    Mach(con) = M;
    
    frac1 = ((gamma+1) * ((M1 * sin(thetaO)).^2) + 2)./((gamma-1) * ((M1 * sin(thetaO)).^2) + 2);
    frac2 = ((gamma+1)/2) + 1./((M1 * sin(thetaO)).^2);
    
    Cp(con) = (thetaO.^2) .* (1 + (frac1 .* log(frac2)));
    Pratio = (2*gamma * M1.^2 * (sin(tau).^2) - (gamma-1))/(gamma+1);
    
else
    % Fundamentals of Aerodynamics, Anderson 2001 & Development of an
    % Aerodynamics Code for the Optimisation of Hypersonic Vehicles Jazra
    % & Smart 2009
    Mn1 = M1*sin(beta);
    
    numer = 1+((gamma-1)/2) * Mn1.^2;
    denom = gamma * (Mn1.^2) - (gamma-1)/2;
    
    Mn2 = (numer./denom).^0.5;
    Mach(con) = (Mn2./(sin(beta-thetaO)));
    Pratio = 1 + ((2*gamma)/(gamma+1)) * ((Mn1.^2) - 1);
    
    % Assumed equation
    Cp(con) = (1/(0.5*gamma*(M1^2))) * (Pratio-1);
end

P(con) = Pratio*Pinf;
method(ID(con)) = 4;

if any(~con)
    
    [Cp(~con),Mach(~con),P(~con)] = newtonian(thetaN,flow);
    method(ID(~con)) = 1;
end