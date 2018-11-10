function [parPos,physicalPos] = conditioning(parPos,cond,varArray,options)

hardstyle = options.hardconditioning;

if hardstyle
    [parPos,physicalPos] = hardcon(parPos,cond,varArray);
else
    [parPos,physicalPos] = versatilecon(parPos,cond,varArray);
end