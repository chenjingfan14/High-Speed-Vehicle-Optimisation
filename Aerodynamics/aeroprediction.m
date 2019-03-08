function [cost,results] = aeroprediction(properties,allPoints,flow,parameters,thetaBetaM,maxThetaBetaM,PrandtlMeyer,options)
%% Aerodynamic prediction

%% Initialisation
% Angle of attack, Mach matrix dimensions and total number of flight states
% to be run
dim = flow.Dim;
runs = flow.Runs;

% Configuration reference parameters
Aref = parameters.Aref;
MAC = parameters.MAC;

Lmat = parameters.L;
Hmat = parameters.H;
Kmat = parameters.LH;

% Takes section properties and outputs point matrices, number of parts,
% whether part is part of body, and prediction methods to be used for every
% section
[points,properties,bodyPart,partType,impactMethod,shadowMethod] = flowfinder(properties);

tot = numel(properties); % Total number of parts

% Initialise total aerodynamic coefficients
initMatrix = zeros(dim);
[Cl,Cd,CN,CA,Cm,L,D,rootMoment,copx] = deal(initMatrix);

wingPressure = cell(dim);

numFoils = 1:sum(~bodyPart); % Total number of aerofoils

% Initialise pressure coffiecient cells for body and aerofoils
[bodyCp,bodyForce,bodyRadLoc,copCell] = deal(cell(dim));
foilCp = cell(dim(1),dim(2),numFoils);

% Option to analyse panels that lay in shadow of prior panels (flow
% direction wise) as shadow instead of impact
shielding = options.Shielding;
viscous = options.Viscous;
control = options.Control;
baseline = options.Baseline;
structures = options.Structure;

for i=1:runs
    %% Outer flight state loop
    
    run = flightstate(flow,i);
    
    alpha = run.alpha;
    Minf = run.Minf;
    delta = run.delta;
    
    Uinf = run.Uinf;
    Pinf = run.Pinf;
    rho = run.rho;
    
    xyAngle = run.planeAngles(1);
    xzAngle = run.planeAngles(2);
    yzAngle = run.planeAngles(3);
        
    if i == 1 || Minf ~= MinfPrev
        % Find maximum compression angle for which attached shock solution exists
        % for free stream Mach number
        loc = Minf == maxThetaBetaM(1,:);
        maxThetaInf = maxThetaBetaM(2,loc);
    end
    
    if viscous && (i == 1 || alpha ~= alphaPrev)
        
%         allPoints = velocitydef(allPoints,run);
    end
    
    if control
        if i == 1
            % Only takes in first non-body properties cell
            points(~bodyPart) = deflectcontrolsufrace(properties{~bodyPart},points(~bodyPart),delta);
        elseif delta ~= deltaPrev
            
            points(~bodyPart) = deflectcontrolsufrace(properties{~bodyPart},points(~bodyPart),delta,deltaPrev);
        end
    end
    
    % Initialise aerodynamic coefficients for each part, every angle of
    % attack
    [partCl,partCd,partCN,partCA,partCm,sumCp] = deal(zeros(tot,1));
    
    xyzCp = zeros(tot,3);
    
    % Used for combining nose/forebody/aftbody flow characteristics
    % together
    bodyOff = 0;
    
    % Rotates points by angle of attack, with mean aerodynamic chord used
    % as origin
    rotPoints = rotate(points,MAC,alpha);
    
    % Initialise shield upper and lower boundaries
    foilCount = 1;
    [yzUpBound,yzLoBound] = deal([0,0]);
    [cellCp,cellForce,cellRadLoc] = deal(cell(1));
    
    count = 1;
    colCp = 0;
    
    %     plotter(rotPoints,"impact")
    
    for j=1:tot
        
        converged = false;
        
        partProp = properties{j};
        part = points(j);
        conical = partProp.Conical;
            
        while ~converged
            %% Inner aerodynamic force and moment coefficient prediction loop
            rotPart = rotPoints(j);
            radialLocation = rotPart.radialLocation;
            
            area = part.area; % Area of each panel
            partdel = rotPart.del;
            
            % Panel unit normals x,y and z
            % Unrotated normals used. Correct?
            unitNorm = part.unitNorm;
            centre = part.centre;
            
            [row,col] = size(area);
            
            %% Characteristics prior to part
            % If current is part of body, give it the previous body part flow
            % characteristics, unless it is the first body part (ie jj = 1).
            % For this and for aerofoils, prior characteristics = freestream
            % conditions
            if bodyPart(j) && j > 1 % Previous body part conditions
                
                cols = (1:col) + bodyOff;
                prevdel = prevBodyAngle(cols);
                prevCp = prevBodyCp(cols);
                prevP = prevBodyP(cols);
                prevMach = prevBodyMach(cols);
                bodyOff = bodyOff + cols(end);
                
            else % Freestream conditions
                
                [prevdel,prevCp] = deal(zeros(1,col));
                prevP = Pinf*ones(1,col);
                prevMach = Minf*ones(1,col);
                
            end
            
            del = [prevdel; partdel]; % Rotated panels inclination to flow
            
            dummy = zeros(row,col);
            Cp = [prevCp; dummy];
            P = [prevP; dummy];
            Mach = [prevMach; dummy];
            
            %% Determine if panel is impacted by flow
            % Centre y,z coordinates of each panel
            
            prev = false(1,col);
            zeroArea = area <= 0;
            logicalFlow = rotPart.flow;
            
            if shielding
                
                cyRot = rotPart.centre(:,Y);
                czRot = rotPart.centre(:,Z);
                
                [impact,shadow] = impactshadow(cyRot,czRot,area,prev,logicalFlow,yzUpBound,yzLoBound);
            else
                impact = [prev; logicalFlow & ~zeroArea];
                shadow = [prev; ~logicalFlow & ~zeroArea];
            end
            
            %% Impact solver
            if any(impact(:))
                
                impactdel = del(impact);
                iMethod = impactMethod(j);
                
                % If using oblique shock/tangent method, initial inclination to
                % flow for attached shock must be less than maximum allowable
                % (assuming front is bluntest part of part for tangent wedge/
                % cone method), otherwise switch to Newtonian Prandtl-Meyer
                if iMethod == 4 && any(abs(partdel(:)) > maxThetaInf)
                    
                    iMethod = 3;
                end
                
                meandel = mean(abs(partdel(1,:)));
                
                if iMethod == 3 && any(meandel > maxThetaInf)
                    
                    iMethod = 2;
                end
                
                meandel = meandel*180/pi;
                
                switch iMethod
                    case 1 % Modified Newtonian
                        
                        [impactCp,impactMach,impactP] = newtonian(partProp,impactdel,meandel,run);
                        
                    case 2 % Modified Newtonian + Prandtl-Meyer expansion
                        
                        [impactCp,impactMach,impactP] = newtonianprandtlmeyer(partProp,del,impact,meandel,Cp,Mach,P,run,PrandtlMeyer,thetaBetaM,maxThetaBetaM,conical);
                        
                    case 3 % Oblique shock + Prandtl-Meyer expansion
                        
                        [impactCp,impactMach,impactP] = obliqueshockprandtl(del,impact,Cp,Mach,P,run,PrandtlMeyer,thetaBetaM,maxThetaBetaM,conical);
                        
                    case 4 % Tangent wedge/cone
                        
                        [impactCp,impactMach,impactP] = tangentobliqueshock(impactdel,run,thetaBetaM,maxThetaBetaM,conical);
                end
                
                Cp(impact) = impactCp;
                Mach(impact) = impactMach;
                P(impact) = impactP;
            end
            
            %% Shadow solver
            if any(shadow(:))
                
                sMethod = shadowMethod(j);
                
                % Initialise shadow panel characteristics as zero arrays
                nrow = sum(shadow(:));
                
                [shadowCp,shadowMach,shadowP] = deal(zeros(nrow,1));
                
                switch sMethod
                    case 1 % Newtonian/High Mach number base pressure
                        % Condition to define whether panel is a base or not,
                        % ie. if angle between it and prior panel is large, or
                        % if it has large shadow inclination
                        
                        %% CHECK THIS
%                         delShadow = del(shadow);
%                         
%                         ddel = delShadow(2:end) - delShadow(1:end-1);
%                         
%                         con = abs(ddel) > (45*pi/180) | abs(delShadow(2:end)) > (80*pi/180);
%                         %%
%                         
%                         shadowCp(con) = - 1/(Minf^2);
%                         shadowP(con) = 0.5*shadowCp(con)*gamma*(Minf^2) + Pinf;
                        
                        % Cp/Mach/P will all be zero (vacuum conditions)
                        % All are initialised to zero thus nothing needs to be
                        % done here
                        
                    case 2 % Prandtl-Meyer expansion
                        [shadowCp,shadowMach,shadowP] = prandtlmeyer(del,shadow,Cp,Mach,P,run,PrandtlMeyer);
                end
                
                Cp(shadow) = shadowCp;
                Mach(shadow) = shadowMach;
                P(shadow) = shadowP;
                
            end
            
            % If panel matrices contain previous conditions, remove to retain
            % only current part characteristics
            [check,~] = size(Cp);
            if check ~= row
                Cp(1,:) = [];
                Mach(1,:) = [];
                P(1,:) = [];
            end
            
            
            %% Structures
            
            if partType(j) == "wing" && structures
               
                wingPressure{i} = P .* -unitNorm;
                [displace, partProp, part] = structure3D(partProp, part, wingPressure{i}, Hmat, Lmat, Kmat);
                
                absDisplace = abs(displace);
                rotPoints(j) = rotate(part,MAC,alpha);
                
%                 if all(absDisplace < 1e-5)
                if all(absDisplace < 1e5)
                    
                    converged = true; 
                else
                    mean(displace)
                    max(absDisplace)
                end
            else
                converged = true;
            end
        end
        
        %% Wing bending moment
        % Only call if part is first aerofoil (wing should always be set up
        % to be first aerofoil in configuration)
        if any(partProp.Name == ["aerofoil","wing"]) && foilCount == 1
            
            rootMoment(i) = wingbending(Cp,rotPoints(j),run);
        end
        
        %% Shielding
        
        if shielding && conical
            
            yzBound = [yzUpBound; yzLoBound];
            [yzUpBound,yzLoBound] = shadowmatrix(cyRot,czRot,yzBound);
        end
        
        %% Characteristics for next part
        % Set final panels as previous characteristics for next part, save
        % Cp for body or aerofoils
        
        if bodyPart(j)
            
            cellCp{count} = Cp;
            cellForce{count} = Cp.*area;
            cellRadLoc{count} = radialLocation;
            
            colCp = colCp + size(cellCp{count},2);
            
            
            if colCp == size(bodyCp{i},2) || j == 1
                
                bodyCp{i} = [bodyCp{i}; cellCp{:}];
                bodyForce{i} = [bodyForce{i}; cellForce{:}];
                bodyRadLoc{i} = [bodyRadLoc{i}; cellRadLoc{:}];
                
                [cellCp,cellRadLoc] = deal(cell(1));
                
                prevBodyAngle = del(end,:);
                prevBodyCp = Cp(end,:);
                prevBodyP = P(end,:);
                prevBodyMach = Mach(end,:);
                
                bodyOff = 0;
                count = 1;
                colCp = 0;
                
            else
                count = count + 1;
            end
        else
            foilCp{i,foilCount} = Cp;
            foilCount = foilCount+1;
        end
        
        %% Calculate total part aerodynamic characteristics
        
        centre = squeeze(reshape(centre,[],1,3));
        unitNorm = reshape(unitNorm,[],1,3);
        
        nx = unitNorm(:,:,1);
        ny = unitNorm(:,:,2);
        nz = unitNorm(:,:,3);
        
        area = reshape(area,[],1);
        Cp = reshape(Cp,[],1);
        
        xyzCp(j,:) = sum(centre.*Cp,1);
        sumCp(j) = sum(Cp);
        
        part.CoP = xyzCp(j)/sumCp(j);
        
        % Part aerodynamic characteristics
        partCl(j) = 2*sum(-((Cp .* area .* nx) * sin(xyAngle)) -...
            ((Cp .* area .* ny) * sin(xzAngle)) -...
            ((Cp .* area .* nz) * sin(yzAngle)))/Aref;
        
        partCd(j) = 2*sum(-((Cp .* area .* nx) * sin(yzAngle)) +...
            ((Cp .* area .* ny) * sin(xzAngle)) -...
            ((Cp .* area .* nz) * sin(xyAngle)))/Aref;
        
        partCm(j) = sum(-(Cp .* area .* nx) +...
            (Cp .* area .* nz))/Aref;
        
        partCN(j) = 2*sum(-((Cp .* area .* nz)))/Aref;
        partCA(j) = 2*sum(-((Cp .* area .* nx)))/Aref;
        
        MinfPrev = Minf;
        alphaPrev = alpha;
        deltaPrev = delta;
    end
    
    %% Friction calculation
    if viscous
%         if partProp.Name == "nose"
%             rot = (alpha*pi/180) + partProp.Rotation;
%         else
%             rot = 0;
%         end

        allPoints = cornervelocities(allPoints, run, Lmat);
        
%         intstreamline(allPoints,run);
        
        % Independent of AoA
        % Cdf = simplefriction(properties,partType,parameters,run);
        Cdf = friction(points,bodyPart,Aref,run);
    else
        Cdf = 0;
    end
    
    %% Total configuration characteristics
    
    % Sum coefficients for flight state, distribute to respective variables
    coeffs = num2cell(sum([partCl,partCd,partCN,partCA,partCm],1));
    
    [Cl(i),Cd(i),CN(i),CA(i),Cm(i)] = deal(coeffs{:});
    
    CA(i) = CA(i) + Cdf;
    Cd(i) = Cd(i) + Cdf;
    
    cop = sum(xyzCp)/sum(sumCp);
    
    copx(i) = cop(1);
    copCell{i} = cop;
    
    L(i) = 0.5 * rho * (Uinf^2) * Cl(i) * Aref;
    D(i) = 0.5 * rho * (Uinf^2) * Cd(i) * Aref;
    %     plotter(rotPoints)
    
end

%% Find averages and save all results
copMaxDiff = max(copx,[],1) - min(copx,[],1);
copMaxDiffbar = mean(copMaxDiff);

Lbar = mean(L(:));
Dbar = mean(D(:));
Clbar = mean(Cl(:));
Cdbar = mean(Cd(:));
Mbar = mean(rootMoment(:));

copBar = copMaxDiffbar/MAC;

results.CN = CN;
results.CA = CA;
results.Cl = Cl;
results.Cd = Cd;
results.Cm = Cm;
results.CoP = copCell;
results.RootMoment = rootMoment;
results.Lift = L;
results.Drag = D;

results.Clbar = Clbar;
results.Cdbar = Cdbar;
results.Mbar = Mbar;
results.Lbar = Lbar;
results.Dbar = Dbar;
results.copBar = copBar;

%% Translate aerodynamic characteristics to cost function values
cost = cost_calculation(results,parameters,options);

results.Cost = cost;