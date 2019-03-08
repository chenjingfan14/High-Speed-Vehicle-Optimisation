function part = cornervelocities(part,flow,L)

% clear all
% close all
%
% flow = flowparameters();
% points = plategen();

Vinf = flow.U;

unitNorm = part.unitNorm;

%% Centre point velocities

T = crossmat(unitNorm, permute(Vinf,[3 1 2]));
Vcentre = crossmat(T, unitNorm);

Vcorner = L * Vcentre(:);

part.FirstPanel = firstPanel;
part.Vc = Vcorner;