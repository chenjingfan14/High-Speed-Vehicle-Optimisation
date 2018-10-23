function [Cp] = inclinedcone(theta,phi,flow)

alpha = flow.alpha;
Minf = flow.Minf;

beta = ((Minf^2)-1)^0.5;

y = log(beta*sin(theta));

x = 0.18145 - 0.20923*y + 0.09092*(y^2) + 0.006876*(y^3)...
    - 0.0062225*(y^4) - 0.000971*(y^5);

Cp0 = (2*(sin(theta)^2))*exp(x);

B = 2.9;

T = sin(B*theta)*tan(theta);

a = [-0.07657, 1.4775, 0.064669;
    0.42339, 0.13241, 0.035871;
    -0.002083, -0.075797, -0.01923;
    0.29898, -0.10011, 0.29589;
    -0.99727, -0.41751, 0.068791;
    -0.039442, 0.10422, 0.063801];

A1 = a(1,1) + a(1,2)*cos(phi) + a(1,3)*cos(2*phi);
A2 = a(2,1) + a(2,2)*cos(phi) + a(2,3)*cos(2*phi);
A3 = a(3,1) + a(3,2)*cos(phi) + a(3,3)*cos(2*phi);
A4 = a(4,1) + a(4,2)*cos(phi) + a(4,3)*cos(2*phi);
A5 = a(5,1) + a(5,2)*cos(phi) + a(5,3)*cos(2*phi);
A6 = a(6,1) + a(6,2)*cos(phi) + a(6,3)*cos(2*phi);

Cpf = Cp0 + ((A1*T) + (A2*T/(Minf^2)) + (A3/(Minf^2)))*(alpha/theta) + ((A4*T)...
    + (A5*T/(Minf^2)) + (A6/(Minf^2)))*((alpha/theta)^2);

Cp = Cpf;