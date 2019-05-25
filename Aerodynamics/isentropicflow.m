function [M2,P,Cp] = isentropicflow(M1,P1,flow)

Minf = flow.Minf;
Pinf = flow.Pinf;
gamma = flow.gamma;

% T2_T1 = (1 + ((gamma-1)/2)*(M1.^2))./(1 + ((gamma-1)/2)*(M2.^2));
% T2 = T2_T1 .* T1;

P2_P1 = T2_T1.^(gamma/(gamma-1));
P2 = P2_P1 .* P1;

M1(id) = M2;

Cp = 2*((P2/Pinf)-1)/(gamma*(Minf^2));