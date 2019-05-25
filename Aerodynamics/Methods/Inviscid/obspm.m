function [Cp,Mach,P] = obliqueshock(del,impact,Cp,Mach,P,prev,flow,maxThetaBetaM)

lb = 0;
gamma = flow.gamma;

delT = [del(1,:) - prev.del; diff(del)];

delT = del;

target = tan(delT);

rows = size(target);
rowArray = 1:rows;

rows = rowArray(any(impact,2));

if any(rows ~= 1)
    
    Mach(1,:) = prev.Mach;
end

for i = rows
    
    impacti = impact(i,:);
    deli = del(i,impacti);
    
    if i == 1
        
        M = prev.Mach';
    else    
        M = Mach(i-1,impacti)';
    end
    
    [Cpi, Mi, Pi] = deal(zeros(size(M)));
    
    ub = halfspace(M,maxThetaBetaM(:,1),maxThetaBetaM(:,3));
    
    tbm = @(B,M,gamma) 2*cot(B) .* ((M.^2) .* (sin(B).^2) - 1)./((M.^2) .* (gamma + cos(2*B)) + 2);
    tbmbi = @(B) tbm(B,M,gamma);
    
    [Beta, ~, flag] = bisection(tbmbi,lb,ub,target(i,:)');
    
    con = isnan(Beta);
    
    if any(con)
        
        [Cpi(con),Mi(con),Pi(con)] = newtonian(deli,flow);
    end
    
    if any(~con)
        
        Mi(~con) = 1/(sin(Beta - delT).^2) .* ((gamma - 1) .* M.^2 * (sin(Beta).^2) + 2)./(2 * gamma * (M.^2) * (sin(Beta).^2) - (gamma - 1));
        Pi(~con) = P(i-1,~con) .* (2 * gamma * (M.^2) .* (sin(Beta).^2) - (gamma - 1))/(gamma + 1);
    end
    
    Cp(i,impacti) = Cpi;
    P(i,impacti) = Pi;
    Mach(i,impacti) = Mi;
    
end