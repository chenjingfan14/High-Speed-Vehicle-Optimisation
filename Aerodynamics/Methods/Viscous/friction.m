function Cdf = friction(partStruct,bodyPart,Aref,flow)

%% Eckert's Reference Temperature method
% Iteration used in Hypersonic and High-Temperature Gas Dynamics

% Constant 0.7 in above ref
Pr = flow.Pr;
R = 287;

% Should Ae values be Ainf??
Pe = flow.Pinf;
Minf = flow.Minf;
Me = flow.Minf;
mu = flow.mu;
a = flow.a;
Ve = Me*a;
Tinf = flow.Tinf;

% Mach 2 X-34
% Tinf = 325 - 272.594;
% Mach 6 X-34
% Tinf = 114 - 272.594;

Te = flow.Tinf;
gamma = flow.gamma;
S = 110;

r = Pr.^(1/3);

Cdf = 0;

for i = 1:numel(partStruct)
    
    part = partStruct(i);
    points = part.Points;
    area = part.area;
    centre = part.centre;
    
    % If part is body part. TODO: Clean this here AND in
    % aeroprediction/flowfinder
    con = bodyPart(i);
    
    % Find leading edge point. If body it must be nose (hence i == 1), or
    % wing/tail. LEx = first row of panels. Else part is a body part that
    % is not the nose, and the LEx is therefore the nose.
    if i == 1 || ~con
        
        LE = (points(1,1:end-1,:) +  points(1,2:end,:))/2;
        
        % If nose, save for future body parts, only need to save one of the
        % xPoints as 
        if con
            
            nose = LE(1,1,:);
        end
    else
        LE = nose;
    end
        
        
    % Distance from leading edge. Must be nose for fore/aftbody
    xyz = centre - LE;
    
    dist = (xyz(:,:,1).^2 + xyz(:,:,2).^2 + xyz(:,:,3).^2).^0.5;
    
    Tt = Te*(1 + (gamma - 1)*((Me^2)/2));
    
    T0 = Tinf*(1 + ((gamma-1)/2)*Minf^2);
    
    Tr = r*(Tt - Te) + Te;
    
    % Pretty sure this is wrong
    % Tw = Tr;
    
    % Assume adiabatic wall T
    Tw = r*(T0 + Te) + Te;
    
    % Original Eckert's reference temperature method
    Tref = Te * (1 + 0.032*(Me^2) + 0.58*((Tw/Te) - 1));
    
    % Thermoelastic Formulation of a Hypersonic Vehicle Control Surface for Control-Oriented Simulation
    Tref = Te + 0.5*(Tw - Te) + 0.22*(Tr - Te);
    
    % Meador-Smart reference temperature method - Laminar
    % Tref = Te * (0.45 + 0.55*(Tw/Te) + 0.16*r*((gamma - 1)/2) * Me^2);
    
    % Meador-Smart reference temperature method - Turbulent
    Tref = Te * (0.5*(1 + (Tw/Te)) + 0.16*r*((gamma - 1)/2) * Me^2);
    
    % Ideal gas law (TRef? P?)
    rhoRef = Pe./(R*Tref);
    % Sutherland's law
    muRef = mu * ((Tref/Te).^1.5) * (Te + S)/(Tref + S);
    
    ReRefx = (rhoRef*Ve*dist)/muRef;
    
    % Accurate up to 10^9 according to above frictionpaper ref
    cf = (0.37./(log(ReRefx).^2.584)) .* area/Aref;
    % Equation from above ref
    % cf = (0.0592./(ReRefx.^0.2)) .* area/Aref;
    % Hypersonic and high-temperature gas dynamics (Meador-Smart)
    cf = 0.02296./((ReRefx).^0.139) .* area/Aref;
    
    cf = 0.02296./((ReRefx).^0.139) .* area/Aref;
    
    % Mangler factor for conical bodies
    if con
        partCdf = 3^0.5 * sum(cf(:));
    else
        partCdf = sum(cf(:));
    end
    
    Cdf = Cdf + partCdf;
    
end