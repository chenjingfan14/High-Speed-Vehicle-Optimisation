function [Cp,Mach,P] = obliqueshock(del,flow,maxThetaBetaM,M1,P1)

lb = 0;
gamma = flow.gamma;

if nargin == 3

    M1 = flow.Minf;
    P1 = flow.Pinf;
end

target = tan(del);

ub = halfspace(M1,maxThetaBetaM(:,1),maxThetaBetaM(:,3));

tbm = @(B,M1,gamma) 2*cot(B) .* ((M1.^2) .* (sin(B).^2) - 1)./((M1.^2) .* (gamma + cos(2*B)) + 2);
tbmbi = @(B) tbm(B,M1,gamma);

[Beta, ~, flag] = bisection(tbmbi,lb,ub,target(:));

con = isnan(Beta);

[Cp,Mach,P] = deal(zeros(size(del)));

if any(con)
    
    [Cp(con),Mach(con),P(con)] = newtonian(del(con),flow);
end

con = ~con;

if any(con)
    
    MachOBS = (1./(sin(Beta - del).^2) .* ((gamma - 1) .* M1.^2 .* (sin(Beta).^2) + 2)./(2 * gamma * (M1.^2) .* (sin(Beta).^2) - (gamma - 1))).^0.5;
    POBS = P1 .* (2 * gamma * (M1.^2) .* (sin(Beta).^2) - (gamma - 1))/(gamma + 1);
    CpOBS = -(P1 - POBS)./(0.5 * gamma * P1 .* M1.^2);
    
    Mach(con) = MachOBS(con);
    P(con) = POBS(con);
    Cp(con) = CpOBS(con);
end