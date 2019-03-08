function [Cp,Mach,P] = tangentobliqueshock(theta,flow,iThetaBetaM,maxiThetaBetaM,conical)

alpha = flow.alpha; Minf = flow.Minf; gamma = flow.gamma; Pinf = flow.Pinf;

dummy = zeros(size(theta));

Cp = dummy;
Mach = dummy;
P = dummy;

[~,dim2] = size(iThetaBetaM);
dim2 = dim2-1;
colArray = 1:dim2;
idCol = Minf == iThetaBetaM(1,:);
idCol = colArray(idCol);

% All panel inclinations to flow for attached shock must be less than
% maximum allowable otherwise switch infringing panels to newtonian flow
con = theta <= maxiThetaBetaM(2,idCol);

thetaO = theta(con);
thetaN = theta(~con);
M = zeros(sum(con),1);

if conical
    % Development of an Aerodynamics Code for the Optimisation of 
    % Hypersonic Vehicles - gives M2 > M1
    % Mix between On hypersonic flow past unyawed cone & Approximate
    % solutions for supersonic flow over wedges and cones

    tau = asin(sin(thetaO) .* (((gamma+1)/2) + (1./((Minf*sin(thetaO)).^2))).^0.5);

%     numer = (Minf^2) .* (cos(tau).^2) .* (1+2*(tau-theta));
%     denom = 1 + ((gamma-1)/2) * (Minf^2) * ((sin(tau).^2) - 2*((tau-theta).^2) .* (cos(theta).^2));
%     
%     M = (numer./denom).^0.5;
    
    % Exact Taylor-Maccoll solution?
    for i = 1:length(tau)
        [~,M(i),~] = solvecone(tau(i),Minf,gamma);
    end
    
    Mach(con) = M;
    
    frac1 = ((gamma+1) * ((Minf * sin(thetaO)).^2) + 2)./((gamma-1) * ((Minf * sin(thetaO)).^2) + 2);
    frac2 = ((gamma+1)/2) + 1./((Minf * sin(thetaO)).^2);
    
    Cp(con) = (thetaO.^2) .* (1 + (frac1 .* log(frac2)));
    Pratio = (2*gamma * Minf.^2 * (sin(tau).^2) - (gamma-1))/(gamma+1);
    
else
    % Fundamentals of Aerodynamics, Anderson 2001 & Development of an
    % Aerodynamics Code for the Optimisation of Hypersonic Vehicles Jazra
    % & Smart 2009
    % 2:end to remove initial Mach number row
    iThetaBetaMvec = iThetaBetaM(2:end,[1 idCol]);

    % Only keep values greater than 0 and below max theta (ie. 
    TBMcon = iThetaBetaMvec(:,2) > 0 & iThetaBetaMvec(:,1) <= maxiThetaBetaM(3,idCol);

    iThetaBetaMvec = iThetaBetaMvec(TBMcon,:);

    absdiff = abs(theta(con)' - iThetaBetaMvec(:,2));
    [~,rows] = min(absdiff,[],1);

    beta = iThetaBetaMvec(rows,1);
    
    Mn1 = Minf*sin(beta);
    
    numer = 1+((gamma-1)/2) * Mn1.^2;
    denom = gamma * (Mn1.^2) - (gamma-1)/2;
    
    Mn2 = (numer./denom).^0.5;
    Mach(con) = (Mn2./(sin(beta-thetaO)));
    Pratio = 1 + ((2*gamma)/(gamma+1)) * ((Mn1.^2) - 1);
    
    % Assumed equation
    Cp(con) = (1/(0.5*gamma*(Minf^2))) * (Pratio-1);
end

P(con) = Pratio*Pinf;

if ~isempty(thetaN)
    [Cp(~con),Mach(~con),P(~con)] = newtonian([],thetaN,0,flow);
end