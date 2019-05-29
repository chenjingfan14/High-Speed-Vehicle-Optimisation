function [Cp,Mach,P,method] = newtonianprandtlmeyer(ID,stream,del,prev,impact,method,Cp,Mach,P,flow,pmFun,maxThetaBetaM)

% Use Newtonian method if inclination is greater than Prandtl-Meyer
% matching angle
matchdel = flow.delq;
newt = impact & del > matchdel;

newt([1 2],:) = impact([1 2],:);

% Newtonian solver
if any(newt(:))
    
    newtdel = del(newt);
    
    [newtCp,newtMach,newtP] = newtonian(newtdel,flow);
    Cp(newt) = newtCp;
    Mach(newt) = newtMach;
    P(newt) = newtP;
    
    method(newt) = 1;
end

% Prandtl-Meyer solver
pm = impact & ~newt;

if any(pm(:))
    
    [Cp,Mach,P] = prandtlmeyer(ID,stream,del,prev,pm,Cp,Mach,P,flow,pmFun,maxThetaBetaM);
    method(pm) = 2;
end

Cp = Cp(impact);
Mach = Mach(impact);
P = P(impact);