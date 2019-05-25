function [Cp,Mach,P] = newtonianprandtlmeyer(ID,stream,del,prev,impact,Cp,Mach,P,flow,pmFun,maxThetaBetaM)

% Use Newtonian method if inclination is greater than Prandtl-Meyer
% matching angle
matchdel = flow.delq;
newt = impact & del > matchdel;

if any(newt(:))
    
    newtdel = del(newt);
    
    [newtCp,newtMach,newtP] = newtonian(partProp,newtdel,meandel,flow);
    Cp(newt) = newtCp;
    Mach(newt) = newtMach;
    P(newt) = newtP;
end

% First panel seeing flow (ie previous is freestream or Cp = 0) must have
% Newtonian/oblique shock solution if it is impacted by flow. Solution to
% panels who's inclination is less than match angle must therefore be found
% by oblique shock
tangFirstPanel = Cp(1,:) == 0 & ~newt(2,:) & impact(2,:);

if any(tangFirstPanel)
    
    tang = false(size(newt));
    tang(2,:) = tangFirstPanel;
    impactdel = del(tang);
    
    [tangCp,tangMach,tangP] = tangentobliqueshock(impactdel,flow,thetaBetaM,maxThetaBetaM);
    Cp(tang) = tangCp;
    Mach(tang) = tangMach;
    P(tang) = tangP;
    
    pm = impact & ~tang & ~newt;
else
    pm = impact & ~newt;
end

if any(pm(:))
    
    [Cp,Mach,P] = prandtlmeyer(ID,stream,del,prev,pm,Cp,Mach,P,flow,pmFun);
end

Cp = Cp(impact);
Mach = Mach(impact);
P = P(impact);