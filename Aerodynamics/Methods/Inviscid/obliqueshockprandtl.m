function [Cp,Mach,P] = obliqueshockprandtl(ID,stream,del,impact,area,prev,Cp,Mach,P,flow,pmFun,maxThetaBetaM)

% NEEDS FIXED

% First panel seeing flow (ie previous is freestream or Cp = 0) must have
% Newtonian/oblique shock solution if it is impacted by flow. Solution to
% panels who's inclination is less than match angle must therefore be found
% by oblique shock

newt = del > flow.delq;
newt = del > halfspace(flow.Minf,maxThetaBetaM(:,1),maxThetaBetaM(:,2));

pm = impact;

if any(newt(:))
    
    [Cp(newt),Mach(newt),P(newt)] = newtonian(del(newt),flow);
    
    pm(newt) = false;
end

count = 1;

% Surface may be extending from a previous one, so if the previous pressure
% coefficient has a non-zero value, shock-expansion must be done from the
% previous panel conditions, not freestream
con = prev.Cp & impact(1,:) & area(1,:) ~= 0 & ~newt(1,:);

id = ID(1,~con)';

while any(con) && count < size(ID,1)
    
    ID(id(con)) = nan;
    id(con) = id(con) + 1;
    
    con = area(id) == 0 | ~impact(id) | newt(id);
    
    count = count + 1;
end

id = id(impact(id));

% Freestream impact (first panel seeing flow)
if any(id)
    
    [Cp(id),Mach(id),P(id)] = obliqueshock(del(id),flow,maxThetaBetaM);
    
    con = ~ismember(ID,id) & ~isnan(ID);
    
    pm(con) = false;
end

% Shock-expansion from previous panel
if any(pm(:))
    
    [Cp(pm),Mach(pm),P(pm)] = prandtlmeyer(ID,stream,del,prev,pm,Cp,Mach,P,flow,pmFun);
end