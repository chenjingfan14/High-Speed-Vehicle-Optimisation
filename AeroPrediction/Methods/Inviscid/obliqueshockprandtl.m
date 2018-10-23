function [Cp,Mach,P] = obliqueshockprandtl(del,impact,Cp,Mach,P,flow,PrandtlMeyer,thetaBetaM,maxThetaBetaM,conical) 

% First panel seeing flow (ie previous is freestream or Cp = 0) must have
% Newtonian/oblique shock solution if it is impacted by flow. Solution to
% panels who's inclination is less than match angle must therefore be found
% by oblique shock

row = 2;
tangFirstPanel = impact(row,:);

if any(tangFirstPanel)
    tang = false(size(impact));
    tang(row,:) = tangFirstPanel;
    impactdel = del(tang);

    [tangCp,tangMach,tangP] = tangentobliqueshock(impactdel,flow,thetaBetaM,maxThetaBetaM,conical);
    Cp(tang) = tangCp;
    Mach(tang) = tangMach;
    P(tang) = tangP;

    pm = impact & ~tang;
else
    pm = impact;
end

if any(pm(:))
    [pmCp,pmMach,pmP] = prandtlmeyer(del,pm,Cp,Mach,P,flow,PrandtlMeyer);
    Cp(pm) = pmCp;
    Mach(pm) = pmMach;
    P(pm) = pmP;
end

Cp = Cp(impact);
Mach = Mach(impact);
P = P(impact);