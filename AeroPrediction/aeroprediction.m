function [cost,aerodynamics] = aeroprediction(properties,flow,parameters,thetaBetaM,maxThetaBetaM,PrandtlMeyer,options)
%% Aerodynamic prediction

%% Initialisation
% Angle of attack, Mach matrix dimensions and total number of flight states
% to be run
dim = flow.dim;
runs = flow.runs;
alphaMat = flow.alpha;
MinfMat = flow.Minf;
% Free stream flow parameters
Pinf = flow.Pinf;
rho = flow.rho;
Uinf = flow.Uinf;

% Configuration reference parameters
Aref = parameters.Aref;
MAC = parameters.MAC;
wingspan = parameters.Wingspan;

% Takes section properties and outputs point matrices, number of parts,
% whether part is part of body, and prediction methods to be used for every
% section
[points,numParts,bodyPart,partType,impactMethod,shadowMethod] = flowfinder(properties);

tot = sum(numParts); % Total number of parts

% Initialise total aerodynamic coefficients
initMatrix = zeros(dim);

[Cl,Cd,CN,CA,Cm,rootMoment,copx] = deal(initMatrix);

numFoils = 1:sum(~bodyPart); % Total number of aerofoils

% Initialise pressure coffiecient cells for body and aerofoils
[bodyCp,bodyForce,bodyRadLoc,copCell] = deal(cell(dim));
foilCp = cell(dim(1),dim(2),numFoils);

% Option to analyse panels that lay in shadow of prior panels (flow
% direction wise) as shadow instead of impact
shielding = options.shielding;
viscous = options.viscous;

for i=1:runs
    %% Outer flight state loop 
    
    run = flightstate(flow,i);
    
    alpha = alphaMat(i);
    Minf = MinfMat(i);
    
    xyAngle = run.planeAngles(1);
    xzAngle = run.planeAngles(2);
    yzAngle = run.planeAngles(3);
    
    if i == 1 || MinfMat(i) ~= MinfMat(i-1)
        % Find maximum compression angle for which attached shock solution exists
        % for free stream Mach number
        loc = Minf == maxThetaBetaM(1,:);
        maxThetaInf = maxThetaBetaM(2,loc);
    end
    
    if viscous && (i == 1 || alphaMat(i) ~= alphaMat(i-1))
        points = velocitydef(points,run);
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
    [partCount,foilCount] = deal(1);
    [yzUpBound,yzLoBound] = deal([0,0]);
    [cellCp,cellForce,cellRadLoc] = deal(cell(1));
    
    count = 1;
    colCp = 0;
    
%     plotter(rotPoints,"impact")
    
    for j=1:tot
        %% Inner aerodynamic force and moment coefficient prediction loop
        partProp = properties{partCount};
        part = points(j);
        rotPart = rotPoints(j);
        conical = partProp.Conical;
        radialLocation = rotPart.radialLocation;
        
        area = part.area; % Area of each panel
        partdel = rotPart.del;
        
        [row,col] = size(area);
        
        Z = (1:col)*3;
        Y = Z - 1;
        X = Y - 1;
        
        % Panel unit normals x,y and z
        % Unrotated normals used. Correct?
        unitNx = part.unitNorm(:,X);
        unitNy = part.unitNorm(:,Y);
        unitNz = part.unitNorm(:,Z);
        
        cx = part.centre(:,X);
        cy = part.centre(:,Y);
        cz = part.centre(:,Z);
        
        %% Characteristics prior to part
        % If current is part of body, give it the previous body part flow
        % characteristics, unless it is the first body part (ie jj = 1).
        % For this and for aerofoils, prior characteristics = freestream
        % conditions
        if bodyPart(partCount) && j > 1 % Previous body part conditions
            
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
            iMethod = impactMethod(partCount);
            
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
            
            sMethod = shadowMethod(partCount);
            
            % Initialise shadow panel characteristics as zero arrays
            nrow = sum(shadow(:));
            
            [shadowCp,shadowMach,shadowP] = deal(zeros(nrow,1));
            
            switch sMethod
                case 1 % Newtonian/High Mach number base pressure
                    % Condition to define whether panel is a base or not,
                    % ie. if angle between it and prior panel is large, or
                    % if it has large shadow inclination
                    
                    %% CHECK THIS
%                     delShadow = del(shadow);
%                     
%                     ddel = delShadow(2:end) - delShadow(1:end-1);
%                     
%                     con = abs(ddel) > (45*pi/180) | abs(delShadow(2:end)) > (80*pi/180);
%                     %%
%                     
%                     shadowCp(con) = - 1/(Minf^2);
%                     shadowP(con) = 0.5*shadowCp(con)*gamma*(Minf^2) + Pinf;
                    
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
        
        %% Wing bending moment
        % Only call if part is first aerofoil (wing should always be set up
        % to be first aerofoil in configuration)
        if part.Name == "aerofoil" && foilCount == 1
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
        
        if bodyPart(partCount)
            
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
        
        cx = reshape(cx,[],1);
        cy = reshape(cy,[],1);
        cz = reshape(cz,[],1);
        unitNx = reshape(unitNx,[],1);
        unitNy = reshape(unitNy,[],1);
        unitNz = reshape(unitNz,[],1);
        area = reshape(area,[],1);
        Cp = reshape(Cp,[],1);
        
        xyzCp(j,:) = sum([cx,cy,cz].*Cp,1);
        sumCp(j) = sum(Cp);
        
        part.CoP = xyzCp(j)/sumCp(j);
        
        % Part aerodynamic characteristics
        partCl(j) = 2*sum(-((Cp.*area.*unitNx)*sin(xyAngle)) - ((Cp.*area.*unitNy)*sin(xzAngle)) - ((Cp.*area.*unitNz)*sin(yzAngle)))/Aref;
        partCd(j) = 2*sum(-((Cp.*area.*unitNx)*sin(yzAngle)) + ((Cp.*area.*unitNy)*sin(xzAngle)) - ((Cp.*area.*unitNz)*sin(xyAngle)))/Aref;
        partCN(j) = 2*sum(-((Cp.*area.*unitNz)))/Aref;
        partCA(j) = 2*sum(-((Cp.*area.*unitNx)))/Aref;
        partCm(j) = sum(-(Cp.*area.*unitNx) + (Cp.*area.*unitNz))/Aref;
        
        if j == sum(numParts(1:partCount))
            partCount = partCount + 1;
        end
        
    end
    
    %% Friction calculation
    if viscous
%         if part.Name == "nose"
%             rot = (alpha*pi/180) + partProp.Rotation;
%         else
%             rot = 0;
%         end
%         rotPoints = cornervelocities(rotPoints,rot,run);
%         
%         intstreamline(rotPoints,run);
        
        % Independent of AoA
        Cdf = simplefriction(properties,partType,parameters,run);
        
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
    
%     plotter(rotPoints,"impact")
    
end

%% Translate aerodynamic characteristics to cost function values
% Usually non-dimensionalised

if any(numFoils)
    L = 0.5*rho*(Uinf.^2).*Cl*Aref;

    copMaxDiff = max(copx,[],1) - min(copx,[],1);
    copMaxDiffbar = mean(copMaxDiff);

    Lbar = mean(L(:));
    Cdbar = mean(Cd(:));
    Mbar = mean(rootMoment(:));

    lift = Lbar/wingspan;
    moment = Mbar/wingspan;
    cop = copMaxDiffbar/MAC;
    
    % Constraint section. Apply penalties if desired values are too
    % high/low
    constrain = [moment,cop];
    minVal = [0,0];
    maxVal = [inf,inf];
    
    penalty = violation(constrain,minVal,maxVal);
    
    cost = [1/lift,Cdbar] + penalty;
    
    % If any cost less than zero, particle swarm will see it as optimal,
    % whereas none of these aerodynamic values should be less than zero
    infCon = cost < 0;
    
    cost(infCon) = inf;
else
    cost = [];
end

aerodynamics.Cl = Cl;
aerodynamics.Cd = Cd;
aerodynamics.Cm = Cm;
aerodynamics.CN = CN;
aerodynamics.CA = CA;
aerodynamics.CoP = copCell;
    
%% To view created configuration, uncomment this
% Leave commented during simulations, otherwise it will plot everything 
% plotter(points,"pause")



















%% Plot Data: Leave commented when running optimisations
% Plot Cl, Cl etc against angle of attack for all flight states

% for i = size(MinfMat,2):-1:1
%     legStr{i} = ['Mach: ' num2str(MinfMat(1,i))];
% end
% 
% figure
% hold on
% grid on
% plot(alphaMat,CN)
% xlabel('Angle of Attack')
% ylabel('C_N')
% legend(legStr)
% hold off
% 
% figure
% hold on
% grid on
% plot(alphaMat,CA)
% xlabel('Angle of Attack')
% ylabel('C_A')
% legend(legStr)
% hold off
% 
% figure
% hold on
% grid on
% plot(alphaMat,Cm)
% xlabel('Angle of Attack')
% ylabel('C_m')
% legend(legStr)
% hold off
% 
% figure
% hold on
% grid on
% plot(alphaMat,Cl)
% xlabel('Angle of Attack')
% ylabel('C_L')
% legend(legStr)
% hold off
% 
% figure
% hold on
% grid on
% plot(alphaMat,Cd)
% xlabel('Angle of Attack')
% ylabel('C_D')
% legend(legStr)
% hold off

% figure
% hold on
% grid on
% plot(alphaMat,rootMoment)
% xlabel('Angle of Attack')
% ylabel('Root Bending Moment')
% legend(legStr)
% hold off

%% Plots Cp at every radial location of body
% Produces hundreds of figures so leave commented during simulations

% data = [AoA', CN, CA, Cm];

% for ii=1:dim
%     
%     meanRadLoc = mean(bodyRadLoc{ii},1);
%     
%     [chordPanels,spanPanels] = size(bodyCp{ii});
%     x = 0:1/(chordPanels-1):1;
%     for jj=1:spanPanels
%         
%         location = round(meanRadLoc(jj),1);
%         
%         figure
%         hold on
%         grid on
%         title(['Body chordwise pressure coefficient at ' num2str(location) '^o'] );
%         plot(x,bodyCp{ii}(:,jj))
%         xlabel('x')
%         ylabel('Cp')
%     end
% end
