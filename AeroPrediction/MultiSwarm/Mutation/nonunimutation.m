function parPos = nonunimutation(parPos,nPop,nVar,varMin,varMax,mutProb,it,maxIt)
%% Non-uniform mutation
% Timestep based mutation. Large mutations possible at beginning of
% simulation, reduces to smaller mutations as PSO iterations increase

con = rand(nPop,nVar) <= mutProb;

a1 = rand(nPop,nVar);
a2 = rand(nPop,nVar);
b = 1;

bool = a1 < 0.5;

id1 = con & bool;
id2 = con & ~bool;

gamma = (a2*(1-(it/maxIt))).^b;

parPos(id1) = parPos(id1) + (varMax(id1)-parPos(id1)).*gamma(id1);
parPos(id2) = parPos(id2) - (parPos(id2)-varMin(id2)).*gamma(id2);