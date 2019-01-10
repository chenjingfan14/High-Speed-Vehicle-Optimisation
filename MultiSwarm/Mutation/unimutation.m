function parPos = unimutation(parPos,nPop,nVar,varMinMat,varMaxMat,mutProb)
%% Uniform mutation
% Creates random number matrix same size as variable matrix. Any random
% number less than mutation probability, the corresponding variable is
% assigned randomly between min and max bounds

con = rand(nPop,nVar) <= mutProb;

parPos(con) = unifrnd(varMinMat(con),varMaxMat(con));