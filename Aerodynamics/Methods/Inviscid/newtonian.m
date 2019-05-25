function [Cp,Mach,P] = newtonian(del,flow)

Minf = flow.Minf;
Pinf = flow.Pinf;
gamma = flow.gamma;
matchMach = flow.Machq;

% if meandel >= 45
%     type = 1;
% else
%     type = 2;
% end

% switch type
%     % Newtonian Impact
%     case 1 % Blunt bodies
        a = 2./(gamma*Minf.^2);
        b = (((gamma+1)^2)*Minf.^2)./((4*gamma*Minf.^2)-(2*(gamma-1)));
        c = gamma/(gamma-1);
        d = (1-gamma+(2*gamma*Minf.^2))./(gamma+1);
        
        K = a.*(((b.^c).*d)-1);
        
%     case 2 % Pointed cones/ogives
%         x = properties.Points(:,:,1);
%         y = properties.Points(:,:,2);
%         z = properties.Points(:,:,3);
%         
%         d = max(max(y(:)),max(z(:)))*2;
%         ln = max(x(:)) - min(x(:));
%         
%         k = 1;
%         
%         a = ((Minf.^2)-1)^0.5;
%         b = sin(atan((0.5^k)*(d/ln)+alpha));
%         
%         K = 2.1+0.5.*((a.*b)^-1);
%         
%     case 3 % Hemisphere
%         a = 1/(gamma*Minf^2);
%         b = ((6*Minf^2)/5)^(gamma/(gamma-1));
%         c = (6/((7*Minf^2)-1))^(1/(gamma-1));
%         
%         K = a*((b*c)-1);
%         
%     case 4 % Real gas flows
%         K = (2*(gamma+1)*(gamma+7))/((gamma+3)^2);
%     case 5 % Aerofoils
%         K = 1.95+(0.3925/((Minf^0.3)*tan(del)));
% end

Cp = K.*sin(del).^2;
P = Pinf.*(1 + (Cp*gamma*Minf.^2)/2);

dim = size(Cp);

Mach = matchMach*ones(dim);