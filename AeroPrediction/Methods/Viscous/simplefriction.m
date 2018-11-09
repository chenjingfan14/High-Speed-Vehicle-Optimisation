function Cdv = simplefriction(properties,partType,parameters,flow)

dim = numel(properties);

M = flow.Minf;
rho = flow.rho;
T = flow.Tinf;
R_air = flow.R;
gamma = flow.gamma;

%% Body Parameters

bodyL = parameters.BodyL; % length of fuselage
bodyW = parameters.BodyW; % radius of the major axis
bodyH = parameters.BodyH; % radius of the minor axis
noseL = parameters.NoseL;
Aref = parameters.Aref;

D_eq = 2 * sqrt(bodyW * bodyH); % equivalant diameter
Sbody = pi * bodyW * bodyH; %frontal area of the elliptical body and nose sections.

fineness = bodyL/D_eq;

u = M * (gamma*R_air*T)^0.5;
q = 0.5*rho*(u^2);

computeBody = true;
for i=dim:-1:1
    
    part = properties{i};
    
    switch partType(i)
        case "nose"
            
        case "forebody"
            
        case {"aftbody","test"}
            if computeBody
                % Calculates friction for FULL body
                Cdw(i) = 3.6/((noseL/D_eq)*(M-1) + 3);
                Cdf(i) = 0.053 * (bodyL/D_eq) * ((M/(q*bodyL))^0.2) * (Sbody/Aref); %note reference area is the total wing planform area.
                computeBody = false;
            end
        case "aerofoil"
            
            partitions = part.Partitions;
            area = part.WetArea;
            chord = part.WetMAC;
            
            for j=partitions:-1:1
                
                % Two wings assumed, need tail condition
                Cdf_wing(j) = 2 * (0.0133/((q*chord(j))^0.2)) * (2*area(j)/Aref);
                
            end
            Cdf(i) = sum(Cdf_wing);
            
    end
    
end

Cdv = sum(Cdf) + sum(Cdw);
