function viewcaller(parPos,cond,varArray,foilData,flow)

[nPop,~] = size(parPos);

% Impose conditions on particles
[partArrays,sectionArray] = partIndexing(cond,varArray);
[~,physicalPos] = conditioning(parPos,cond,varArray);

% Assign 2D section matrices to particles. Foils variable = section indices
sections = foilData(physicalPos(:,sectionArray));
for i=1:nPop
    [assemblyProperties,~,MAC,flag] = particlecreator(physicalPos(i,:),partArrays,sections(i,:));
    if ~flag
        particleviewer(assemblyProperties,MAC,flow,i);
    end
    
    if rem(i,50) == 0
        i
    end
    
end