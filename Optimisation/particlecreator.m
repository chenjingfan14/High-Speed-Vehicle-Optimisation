function [cost,results,assemblyProperties,fricData,parPos,parameters] = particlecreator(parPos,configPos,varArray,sections,options)
%% Assembles arbitrary configuration
% Will fuse together multiple aerofoil sections with body. Will also work
% for aerofoil alone or body alone configurations

Wing = options.Wing;
Aft = options.Aft;
Control = options.Control;
Structure = options.Structure;
Viscous = options.Viscous;

%% Create body
% Check if body is part of configuration, if so build it, if not set body
% parts to empty. Currently aftbody defines body existence, and if it does
% nose and forebody must also exist

if Aft
    
    geom = ["UpperLength","yUpperRad","yBotRatio","zUpperRad",...
        "SideLength","zLowerRad"];
    
    % Assign correct variables to their definitions
    noseLength = configPos(varArray == "NoseLength");
    noseRad = configPos(varArray == "NoseRad");
    zNoseOff = configPos(varArray == "zNoseOffset");
    foreBodyLength = configPos(varArray == "ForeLength");
    aftbodyGeom = configPos(any(varArray == geom,2));
    aftBodyLength = configPos(varArray == "AftLength");
    
    noseForeLength = noseLength + foreBodyLength;
    
    % Create aft body portion of body
    aftbody = arbitraryfuse([aftbodyGeom,aftBodyLength],noseForeLength);
    
    [~,dim_f,~] = size(aftbody.Points);
    radial = squeeze(aftbody.Points(1,:,:));
    
    aftbodyHeight = aftbody.Height;
    aftbodyWidth = aftbody.Width;

    parameters.ForeL = noseForeLength;
    parameters.AftL = aftBodyLength;
    parameters.BodyL = noseForeLength + aftBodyLength;
    parameters.BodyW = aftbodyHeight;
    parameters.BodyH = aftbodyWidth;
    
    % Create nose and link together with aftbody by interpolating to
    % create forebody
    
    % If zero length or radius, nose must be pointed, therefore empty as
    % first forebody points will equal singular nose point
    if ~all([noseLength,noseRad])
        
        Nose = '';
        nosePoints = zeros(1,dim_f,3);
    else
        Nose = nose([noseRad,noseLength,zNoseOff],dim_f,noseForeLength);
        nosePoints = Nose.Points;
    end
    
    forebody = foreGen(foreBodyLength,nosePoints,radial,0.5);
    
    body = [nosePoints(1:end-1,:,:); forebody.Points(1:end-1,:,:); aftbody.Points];
    bodyPoints = normals(body,Viscous,options.Quad);
    
    if Viscous
        
        bodyPoints.L = fsiinterpolation(options.RBFfun,bodyPoints);
    
        fricData{1} = bodyPoints;
    else
        fricData{1} = [];
    end
    
else % Set all body parts to empty
    
    aftbodyHeight = 0;
    aftbodyWidth = 0;
    
    aftbody = ''; forebody = ''; Nose = '';
end

%% Create aerofoils
% If no aerofoils in configuration, loop will be skipped as dim = 0
% Count back from number of aerofoils for preallocation purposes

% Update later to include tail
if Wing
    
    wingDim = 1;
    chordDisc = options.ChordDisc;
end

for i=wingDim:-1:1
    
    dihedralVar = varArray == "Dihedral";
    chordVar = varArray == "Chord";
    TEsweepVar = varArray == "TESweep";
    semispanVar = varArray == "Semispan";
    xOffsetVar = varArray == "xOffset";
    zOffsetVar = varArray == "zOffset";
    
    % Assign variables to their physical attribute
    dihedral = configPos(dihedralVar);
    chord = configPos(chordVar);
    TEsweep = configPos(TEsweepVar);
    semispan = configPos(semispanVar);
    xOffset = configPos(xOffsetVar);
    zOffset = configPos(zOffsetVar);
    
    % As well as their "particle" attribute which may need to be updated if
    % configuration cannot be created
    parDihedral = parPos(dihedralVar);
    parChord = parPos(chordVar);
    parTESweep = parPos(TEsweepVar);
    parSemispan = parPos(semispanVar);
    parxOffset = parPos(xOffsetVar);
    parzOffset = parPos(zOffsetVar);
    
    % Ensure initial semispan is larger than maximum body radius
    minFirstSemispan = max([aftbodyHeight,aftbodyWidth]/2);
    
    % Stop semispan getting too large, shouldn't need to be larger than
    % maximum body radius
    maxFirstSemispan = max(aftbodyHeight,aftbodyWidth)*2;
    
    if semispan(1) < minFirstSemispan
        
        [semispan(1),parSemispan(1)] = deal(minFirstSemispan);
    end
    
    % Create aerofoil
    if Control
%         control = configPos(varArray == ");
        liftSurface = wingtail(dihedral,semispan,chord,TEsweep,sections,chordDisc,control);
    else 
        liftSurface = wingtail(dihedral,semispan,chord,TEsweep,sections,chordDisc);
    end
        
    if Aft
        
        % Factor to determine how much the wing/tail is shifted in 
        % necessary direction
        fixFactor = 0.1;
        
        % Number of historical failures to take into account
        nHist = 5;
        
        % Condition to keep shift factor from reducing every failed attempt
        stayOut = -nHist;
        
        % Initialise historical failure reason matrix
        reasonHistory = zeros(nHist,1);
        
        % First attempt at joining wing/tail and body
        attempt = 1;
        
        % Join body and aerofoil together
        [liftSurface,aftbody,success,reason] = bodyfoil(aftbody,liftSurface,[xOffset,zOffset]);
        
        %% If unsuccessful, while loop to fix until it is successful
        % Reason 1: Wing fore of body (increase x offset)
        % Reason 2: Wing aft of body (decrease x offset)
        % Reason 3: Wing below body (increase z offset)
        % Reason 4: Wing above body (decrease z offset)
        % Reason 5: First wing partition within body (increase semispan)
        
        while ~success
            
            % Reason config failed to create
            reasonHistory(1) = reason;
            
            % If while loop is bouncing between fail reasons
            bouncing = all(diff(reasonHistory) ~= 0);
            
            % If history of failure reasons is filled, and loop is bouncing
            % between reasons, and loop hasn't reduced fixFactor recently: 
            % reduce fixFactor & reset stayOut
            if attempt > nHist && bouncing && stayOut >= 0
                
                fixFactor = fixFactor/2;
                stayOut = -nHist;
            end
            
            if isempty(reason)
                % Use to debug if configuration cannot be assembled
                % reason
            
                
            % Fixing for specific reasons. Will always go inside here if
            % wing/tail is before aftbody, or semispan must be increased.
            % Other failure reasons can be reduced by last resort below
            elseif reason == 1 || (reason == 5 && semispan(1) < maxFirstSemispan) || attempt < 10
                switch reason
                    case {1,2}
                        
                        % Choose largest to define wing/tail shift
                        % Chord can't be negative which keeps things nice
                        % here (in terms of signs)
                        physical = [xOffset,chord(1)];
                        particle = [parxOffset,parChord(1)];
                        
                        [~,use] = max(physical);
                        
                        if reason == 1 % Increase offset
                            xOffset = xOffset + fixFactor*physical(use);
                            parxOffset = parxOffset + fixFactor*particle(use);
                        
                        else % Reduce offset
                            xOffset = xOffset - fixFactor*physical(use);
                            parxOffset = parxOffset - fixFactor*particle(use);
                        end
                        
                    case {3,4}
                        
                        % Should not come in here for vertical tails
                        
                        % zOffset defined from centre of aftbody, so
                        % negative will bring wing towards centre
                        
                        zOffset = zOffset - fixFactor*zOffset;
                        parzOffset = parzOffset - fixFactor*parzOffset;
                        
                        if zOffset == 0 || (attempt > 5 && reason == 4)
                            zOffset = zOffset + fixFactor;
                            parzOffset = parzOffset + fixFactor;
                        end
                        
                    case 5 % Increase semispan
                        
                        semispan(1) = semispan(1) + 0.1 * semispan(1);
                        parSemispan(1) = parSemispan(1) + 0.1 * parSemispan(1);
                            
                        if Control
                            liftSurface = wingtail(dihedral,semispan,chord,TEsweep,sections,chordDisc,control);
                        else
                            liftSurface = wingtail(dihedral,semispan,chord,TEsweep,sections,chordDisc);
                        end
                end
                
            else
                
                % If wing cannot be combined with body, quit and give inf
                % values to cost functions
                cost = inf;
                assemblyProperties = [];
                parameters = [];
                return
            end
            
            % Re-try merge
            [liftSurface,aftbody,success,reason] = bodyfoil(aftbody,liftSurface,[xOffset,zOffset]);
            
            % Increase number of attempts and stayOut condition
            attempt = attempt + 1;
            stayOut = stayOut + 1;
            
            % Shift current reason down one in the history array to make 
            % room for new reason
            reasonHistory = circshift(reasonHistory,1);
        end
        
        %%
        liftSurface.Offset = [xOffset,zOffset];
        
        % No need to feedback if successful first time
        if attempt > 1
            % Feedback particle and physical positions to optimiser
            parPos(dihedralVar) = parDihedral;
            parPos(semispanVar) = parSemispan;
            parPos(chordVar) = parChord;
            parPos(TEsweepVar) = parTESweep;
            parPos(xOffsetVar) = parxOffset;
            parPos(zOffsetVar) = parzOffset;
            
            configPos(dihedralVar) = dihedral;
            configPos(semispanVar) = semispan;
            configPos(chordVar) = chord;
            configPos(TEsweepVar) = TEsweep;
            configPos(xOffsetVar) = xOffset;
            configPos(zOffsetVar) = zOffset;
        end
    end
    
    if Structure && i == 1
        
        wingbox.Location = [0.2,0.8]';
        wingbox.SparThick = 0.05;
        wingbox.SkinThick = 0.02;
        
        liftSurface = wingstructure(liftSurface,wingbox);
        
    end
    
    liftSurfPoints = normals(liftSurface.Points,Viscous,options.Quad);
    
    if Structure
    
        [L, H, K] = fsiinterpolation(options.RBFfun,liftSurfPoints,liftSurface.Skin);
    
        liftSurface.L = L;
        liftSurface.H = H;
        liftSurface.K = K;

        liftSurfPoints.L = L;
    
    elseif Viscous
       
        L = fsiinterpolation(options.RBFfun,liftSurfPoints);
        liftSurfPoints.L = L;
    end
    
    if Viscous
        
%         j = 1;
        j = length(fricData) + 1;
        fricData{j,1} = liftSurfPoints;
    else
        fricData = [];
    end
    
    liftSurfaces(i) = liftSurface;    
end

if Wing
    
    wing = liftSurfaces(1);
    
    parameters.Chord = chord;
    parameters.Semispan = semispan;
    parameters.Dihedral = dihedral;
    
    if Control
        parameters.Control = control;
    end
    
    sumArea = sum(wing.Area);
    
    parameters.MAC = sum(wing.WetArea.*wing.WetMAC)/sumArea;
    parameters.Aref = sumArea*2;
    parameters.Wingspan = sum(wing.Span);
    parameters.Sweep = wing.LESweep;
    parameters.ThicknessND = wing.ThicknessND; 
    parameters.Thickness = wing.Thickness;
    
else
    
    liftSurfaces = '';
    
    parameters.Aref = aftbody.Area;
    parameters.MAC = aftbody.Length/2;
    parameters.BodyW = aftbody.Width;
    parameters.wingspan = [];
    
end

% X-34 parameters
% parameters.Aref = 33.213;
% parameters.MAC = 4.4323;

%% Put assembly into single cell, delete any empty parts
assemblyProperties = {Nose,forebody,aftbody,liftSurfaces};
assemblyProperties(strcmp(assemblyProperties,'')) = [];

[cost,results] = aeroprediction(assemblyProperties,fricData,parameters,options);