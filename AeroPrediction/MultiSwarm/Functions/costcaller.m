function [parFit,count] = costcaller(costFun,nPop,nFun,parPos,partArrays,sections,flow,thetaBetaM,maxThetaBetaM,PrandtlMeyer,options)

% Initialise fitness matrix
parFit = zeros([nPop, nFun]);

count = 0; % Count successful configuration creations
if options.parallel
    parfor i=1:nPop
        % Attempt to create configuration
        [parProperties,parameters,flag] = particlecreator(parPos(i,:),partArrays,sections(i,:));
        if flag % If unsuccessful, set fitness values to infinity
            parFit(i,:) = ones(1,nFun)*inf;
        else % Else use cost function to calculate fitness
            parFit(i,:) = feval(costFun,parProperties,flow,parameters,thetaBetaM,maxThetaBetaM,PrandtlMeyer,options);
            count = count+1;
        end
    end
else % Same as above but for single-core processing
    for i=1:nPop
        [parProperties,parameters,flag] = particlecreator(parPos(i,:),partArrays,sections(i,:));
        if flag
            parFit(i,:) = ones(1,nFun)*inf;
        else
            parFit(i,:) = feval(costFun,parProperties,flow,parameters,thetaBetaM,maxThetaBetaM,PrandtlMeyer,options);
            count = count+1;
        end
    end
end

% For any particle with cost function values with imaginary parts or not a 
% number, set all particle cost function values to infinity
con = any(~isreal(parFit) & imag(parFit) > 0 | isnan(parFit),2);
parFit(con,:) = inf;