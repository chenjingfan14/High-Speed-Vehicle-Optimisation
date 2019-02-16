function [parFit,parPos] = costcaller(costFun,nPop,nFun,parPos,cond,varArray,n,foilData,flow,thetaBetaM,maxThetaBetaM,PrandtlMeyer,options)

% Initialise fitness matrix
parFit = zeros([nPop, nFun]);

% Impose conditions on particles
Bezier = options.Bezier;

[parPos,physicalPos] = contrans(parPos,nPop,cond,varArray,options);

if Bezier
    % Create 2D aerofoil section Bezier curves from control points
    sectionArray = varArray == "Bezier";
    sectionPos = physicalPos(:,sectionArray);
    
    sections = Bezier3(sectionPos,n,foilData,nPop);
else
    % Assign 2D section matrices to particles. Foils variable = section indices
    sectionArray = varArray == "Section";
    sectionPos = physicalPos(:,sectionArray);
    
    zero = sectionPos == 0;
    sectionPos(zero) = 1;
    sections = foilData(sectionPos);
end

if options.Parallel
    % options.parallel
    parfor i=1:nPop
        
        % Attempt to create configuration
        [parProperties,parPos(i,:),parameters,flag] = particlecreator(parPos(i,:),physicalPos(i,:),varArray,sections(i,:),options);
        
        if flag
            
            parFit(i,:) = inf;
        else
        % End use cost function to calculate fitness
        parFit(i,:) = feval(costFun,parProperties,flow,parameters,thetaBetaM,maxThetaBetaM,PrandtlMeyer,options);
        end
    end
else % Same as above but for single-core processing
    for i=1:nPop
        
        % Attempt to create configuration
        [parProperties,parPos(i,:),parameters,flag] = particlecreator(parPos(i,:),physicalPos(i,:),varArray,sections(i,:),options);
                
        if flag
            
            parFit(i,:) = inf;
        else
        % End use cost function to calculate fitness
        parFit(i,:) = feval(costFun,parProperties,flow,parameters,thetaBetaM,maxThetaBetaM,PrandtlMeyer,options);
        end
    end
end

% For any particle with cost function values with imaginary parts or not a 
% number, set all particle cost function values to infinity
con = any(~isreal(parFit) & imag(parFit) > 0 | isnan(parFit),2);
parFit(con,:) = inf;