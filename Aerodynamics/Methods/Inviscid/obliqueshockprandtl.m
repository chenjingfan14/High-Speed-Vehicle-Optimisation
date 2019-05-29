function [Cp,Mach,P,method] = obliqueshockprandtl(ID,stream,del,impact,method,prev,Cp,Mach,P,flow,pmFun,maxThetaBetaM)

% NEEDS FIXED

% First panel seeing flow (ie previous is freestream or Cp = 0) must have
% Newtonian/oblique shock solution if it is impacted by flow. Solution to
% panels who's inclination is less than match angle must therefore be found
% by oblique shock
newt = del > halfspace(flow.Minf,maxThetaBetaM(:,1),maxThetaBetaM(:,2));

pm = impact;

if any(newt(:))
    
    [Cp(newt),Mach(newt),P(newt)] = newtonian(del(newt),flow);
    
    pm(newt) = false;
    method(newt) = 1;
end

% Surface may be extending from a previous one, so if the previous pressure
% coefficient has a non-zero value, shock-expansion must be done from the
% previous panel conditions, not freestream
shock = false(size(ID));

shock(1,:) = prev.Cp & impact(1,:) & ~newt(1,:);
con = impact(1,:) == 0;
shock(2,con) = prev.Cp(con) & impact(2,con) & ~newt(2,con);

id = ID(shock);

% Freestream impact (first panel seeing flow)
if any(id)
    
    [Cp(id),Mach(id),P(id),method] = obliqueshock(del(id),ID(id),method,flow,maxThetaBetaM);
    
    pm(id) = false;
end

% Shock-expansion from previous panel
if any(pm(:))
    
    [Cp,Mach,P] = prandtlmeyer(ID,stream,del,prev,pm,Cp,Mach,P,flow,pmFun,maxThetaBetaM);
    
    method(pm) = 2;
end