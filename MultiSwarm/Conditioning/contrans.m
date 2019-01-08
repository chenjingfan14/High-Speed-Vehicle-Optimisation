function [parPos,physicalPos] = contrans(parPos,nPop,cond,varArray,options)

hardstyle = options.Hardtransform;

parPos = conditioning(parPos,nPop,cond);

if hardstyle
    physicalPos = hardtransform(parPos,cond,varArray);
else
    physicalPos = versatiletransform(parPos,cond,varArray);
end