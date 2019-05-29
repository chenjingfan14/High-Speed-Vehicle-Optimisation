function [Cp,Mach,P] = prandtlmeyer(ID,stream,del,prev,pm,Cp,Mach,P,flow,pmFun,maxThetaBetaM)

ID(~pm) = 0;
trace = ismember(stream,ID);

Minf = flow.Minf;
Pinf = flow.Pinf;
gamma = flow.gamma;
vmax = (pi/2)*(((gamma+1)/(gamma-1))^0.5 - 1);

[row,col] = size(stream);

dTheta = [del(1,:); diff(del)];
dTheta = dTheta(stream);

rowCon = any(trace,2);
rowArray = 1:row;
rowi = rowArray(rowCon);
colArray = (1:col)';

Cpbase = basepressure(Minf,gamma);

for i = rowi
    
    j = trace(i,:);
    colj = colArray(j);
    
%     dThetai = dTheta(i-1,j)';
    dThetai = dTheta(i,j)';
    
    id = i + (colj-1)*row;
    
    if i == 1 || ~any(Mach(id-1))
        
        P1 = prev.P(j)';
        M1 = prev.Mach(j)';
    else
        P1 = P(id-1);
        M1 = Mach(id-1);
        
        %% Should be shockwaves here
        % Although current vehicle gen should not make configs that have
        % contracting surfaces, may just be zero area panels causing negatives
    end
    
    if any(M1 < 1)
        
        % error('Subsonic flow (M1 < 1)')
    end
    
    % Any undefined Mach/pressures set to freestream
    M1(M1 == 0) = Minf;
    P1(P1 == 0) = Pinf;
    
    % Absolute here for first panels which may have negative inclination to
    % the flow. Otherwise should always be positive?
    vMp1 = pmFun(M1,gamma);
    vMp2 = vMp1 + abs(dThetai);
    
    pmFun2 = @(M1) pmFun(M1,gamma);
    
    [M2,~,flag] = bisection(pmFun2,M1,1000,vMp2);
    
    if any(M2 < M1)
        
        error('Prandtl-Meyer exception (M2 < M1)')
    end
    
    shock = dThetai > 1e-10;
    
    if any(shock)

        [CpShock,MachShock,PShock] = obliqueshock(dThetai(shock),ID(id(shock)),[],flow,maxThetaBetaM,M1(shock),P1(shock));
        
        Mach(id(shock)) = MachShock;
        P(id(shock)) = PShock;
        Cp(id(shock)) = CpShock;
    end
    
    base = isnan(M2);
    
    if any(base)
        
        idbase = id(base);
        
        Mach(idbase) = Minf;
        P(idbase) = Pinf;
        Cp(idbase) = Cpbase;
        
%         Mach(idbase) = 0; 
%         P(idbase) = 0;
%         Cp(idbase) = 0;
    end

    pmi = ~base & ~shock;
    
    if any(pmi)
        
        id = id(pmi);
   
        T2_T1 = (1 + ((gamma-1)/2)*(M1.^2))./(1 + ((gamma-1)/2)*(M2.^2));
    %     T2 = T2_T1 .* T1;

        P2_P1 = T2_T1.^(gamma/(gamma-1));
        P2 = P2_P1 .* P1;

        Cp2 = 2 * ((P2/Pinf) - 1)/(gamma * (Minf^2));

        Mach(id) = M2(pmi);
        P(id) = P2(pmi);
        Cp(id) = Cp2(pmi);
    end
    
    %% Keep? Destroy?
    con = del(id) > flow.delq;
    
    if any(con)
        
        [Cp(id(con)),Mach(id(con)),P(id(con))] = newtonian(del(id(con)),flow);
    end
    
    %% Alterations based on max turning angle
    % If turning angle greater than max, use different method to find panel
    % characteristics
%     con = dThetai > vmax - vMp1;
%     
%     if any(con)
%         %     % Base pressure
%         %     Cp(i,con) = 1/(Minf^2);
%         %     % Newtonian
%         %     Cp(i,con) = 0;
%         %     Mach(i,con) = Minf;
%         %     P(i,con) = 0.5*Cp(i,con)*gamma*(Minf^2) + Pinf;
% 
%         % Vacuum
%         Cp(i,con) = 0;
%         Mach(i,con) = 0;
%         P(i,con) = 0;
%     end
end

% Cp = Cp(pm);
% Mach = Mach(pm);
% P = P(pm);