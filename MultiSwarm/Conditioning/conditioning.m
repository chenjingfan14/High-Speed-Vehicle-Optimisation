function [parPos,physicalPos] = conditioning(parPos,nPop,cond,varArray,options)

hardstyle = options.hardconditioning;

if hardstyle
    [parPos,physicalPos] = hardcon(parPos,nPop,cond,varArray);
else
    [parPos,physicalPos] = versatilecon(parPos,nPop,cond,varArray);
end