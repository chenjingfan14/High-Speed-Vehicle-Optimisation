function parPos = unimutation(parPos,nPop,nVar,varMin,varMax,mutProb)
%% Uniform mutation
% Creates random number matrix same size as variable matrix. Any random
% number less than mutation probability, the corresponding variable is
% assigned randomly between min and max bounds

[dim,~] = size(parPos);

varMinMat = repmat(varMin,dim,1);
varMaxMat = repmat(varMax,dim,1);

con = rand(nPop,nVar) <= mutProb;

parPos(con) = unifrnd(varMinMat(con),varMaxMat(con));