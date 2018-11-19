function viewcaller(parPos,cond,varArray,foilData,n,flow,options)

[nPop,~] = size(parPos);

% Impose conditions on particles
[partArrays,sectionArray] = partIndexing(cond,varArray);
[~,physicalPos] = conditioning(parPos,cond,varArray,options);

% Assign 2D section matrices to particles. Foils variable = section indices
if options.Bezier
    sections = Bezier3(physicalPos(:,sectionArray),n,foilData,nPop);
else
    sections = foilData(physicalPos(:,sectionArray));
end

for i=1:nPop
    success = false;
    attempt = 1;
    
    % Attempt to create configuration
    while ~success
        [parProperties,parPos(i,:),physicalPos(i,:),~,success] = particlecreator(parPos(i,:),physicalPos(i,:),partArrays,sections(i,:),attempt);
        attempt = attempt + 1;
    end
    
    particleviewer(parProperties,i);
    
    if rem(i,50) == 0
        i
    end
    
end