function [thetaBetaM,maxThetaBetaM] = thetabetamachcurves(Mach,beta)

if nargin < 1
    
    Mach = 1.01:0.01:20;
end

if nargin < 2
    
    beta = 0:0.0001:pi/2;
end

dim = length(Mach);
thetaMax = zeros(size(Mach));

id = zeros(size(Mach));
theta = zeros(dim,length(beta));

for i = 1:length(Mach)

    M = Mach(i);

    gamma = 1.4;

    theta(i,:) = atan(2*cot(beta).*((M^2*(sin(beta).^2)-1)./(M^2*(gamma+cos(2*beta))+2)));

    [thetaMax(i),id(i)] = max(theta(i,:));

end

betaMax = beta(id);
beta = [NaN; beta'];
thetaM = [Mach; theta']; 
thetaBetaM = [beta, thetaM];

maxThetaBetaM = [Mach; thetaMax; betaMax]';

save('thetaBetaCurves','thetaBetaM','maxThetaBetaM')