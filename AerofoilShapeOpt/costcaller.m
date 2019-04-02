function [parFit,parPos] = costcaller(costFun,nPop,nFun,parPos,cond,varArray,options)

% Initialise fitness matrix
parFit = zeros([nPop, nFun]) + inf;

% Impose conditions on particles
aerofoilMethod = options.AerofoilMethod;
n = options.WingPartitions;

% Transformation/conditioning might not always be necessary for different
% optimisations
    
[parPos,physicalPos] = contrans(parPos,nPop,cond,varArray,options);

if aerofoilMethod == "BP3434"

    uDisc = options.UpperDisc;
    lDisc = options.LowerDisc;
    
    sections = BP3434(parPos,varArray,n,nPop,lDisc,uDisc);
    
elseif aerofoilMethod == "Bezier"
    
    uDisc = options.UpperDisc;
    lDisc = options.LowerDisc;
    
    % Create 2D aerofoil section Bezier curves from control points
    sectionArray = varArray == "Bezier";
    sectionPos = physicalPos(:,sectionArray);
    
    sections = Bezier3(sectionPos,n,options.BezierControlPoints,nPop,lDisc,uDisc);
    
elseif aerofoilMethod == "BezierTC"
    
    uDisc = options.UpperDisc;
    lDisc = options.LowerDisc;
    
    % Create 2D aerofoil section Bezier curves from control points
    sectionArray = varArray == "Bezier";
    sectionPos = physicalPos(:,sectionArray);
    
    sections = BezierTC(sectionPos,n,options.BezierControlPoints,nPop,lDisc,uDisc);
    
elseif aerofoilMethod == "Preloaded"
    
    % Assign 2D section matrices to particles. Foils variable = section indices
    sectionArray = varArray == "Section";
    sectionPos = physicalPos(:,sectionArray);
    
    zero = sectionPos == 0;
    sectionPos(zero) = 1;
    sections = options.foilData(sectionPos);
end

if options.Parallel
    
    parfor i=1:nPop

        % Use cost function to calculate fitness
        % Aerostructural cost call
%         parFit(i,:) = feval(costFun,parPos(i,:),physicalPos(i,:),varArray,sections(i,:),options);
        % Aerofoil shape comparison cost call
        parFit(i,:) = feval(costFun,sections{i,:},options.Compare);
    end
else % Same as above but for single-core processing
    for i=1:nPop
        
%         parFit(i,:) = feval(costFun,parPos(i,:),physicalPos(i,:),varArray,sections(i,:),options);
        if ~isempty(sections{i})
            
            parFit(i,:) = feval(costFun,sections{i,:},options.Compare);
        end
    end
end

% For any particle with cost function values with imaginary parts or not a 
% number, set all particle cost function values to infinity
con = any(~isreal(parFit) & imag(parFit) > 0 | isnan(parFit),2);
parFit(con,:) = inf;