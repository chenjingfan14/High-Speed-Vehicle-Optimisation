function [assemblyProperties,parPos,parameters] = particlecreator(parPos,configPos,varArray,sections,options)
%% Assembles arbitrary configuration
% Will fuse together multiple aerofoil sections with body. Will also work
% for aerofoil alone or body alone configurations

Wing = options.Wing;
Aft = options.Aft;
Fore = options.Fore;
Nose = options.Nose;
Control = options.Control;
Structure = options.Structure;

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
    
    parameters.ForeL = noseForeLength;
    parameters.AftL = aftBodyLength;
    parameters.BodyL = noseForeLength + aftBodyLength;
    parameters.BodyW = aftbody.Width;
    parameters.BodyH = aftbody.Height;
    
else % Set all body parts to empty
    aftbody = ''; forebody = ''; Nose = '';
end

%% Create aerofoils
% If no aerofoils in configuration, loop will be skipped as dim = 0
% Count back from number of aerofoils for preallocation purposes

% Update later to include tail
if Wing
    wingDim = 1;
end

for i=wingDim:-1:1
    
    dihedralVar = varArray == "Dihedral";
    chordVar = varArray == "Chord";
    sweepVar = varArray == "LESweep";
    semispanVar = varArray == "Semispan";
    xOffsetVar = varArray == "xOffset";
    zOffsetVar = varArray == "zOffset";
    
    % Assign variables to their physical attribute
    dihedral = configPos(dihedralVar);
    chord = configPos(chordVar);
    sweep = configPos(sweepVar);
    semispan = configPos(semispanVar);
    xOffset = configPos(xOffsetVar);
    zOffset = configPos(zOffsetVar);
    
    % As well as their "particle" attribute which may need to be updated if
    % configuration cannot be created
    parDihedral = parPos(dihedralVar);
    parChord = parPos(chordVar);
    parSweep = parPos(sweepVar);
    parSemispan = parPos(semispanVar);
    parxOffset = parPos(xOffsetVar);
    parzOffset = parPos(zOffsetVar);
    
    % If first semispan = 0, give starting value so that wing can be
    % created
    if semispan(1) == 0
        semispan(1) = 0.1;
        parSemispan(1) = 0.1;
    end
    
    % Create aerofoil
    if Control
%         control = configPos(varArray == ");
        liftSurface(i) = wingtail(dihedral,semispan,chord,sweep,sections,control);
    else 
        liftSurface(i) = wingtail(dihedral,semispan,chord,sweep,sections);
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
        
        % Stop semispan getting too large, shouldn't need to be larger than
        % maximum body radius
        maxFirstSemispan = max(aftbody.Height,aftbody.Width)*2;
        
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
            % reduce fixFactor & reset stayOut EDIT: Trying to reduce chord
            % instead
            if attempt > nHist && bouncing && stayOut >= 0
                
%                 chord(1) = chord(1) - 0.1*chord(1);
%                 parChord(1) = parChord(1) - 0.1*parChord(1);
%                 
%                 if Control
%                     liftSurface(i) = wingtail(dihedral,semispan,chord,sweep,sections,control);
%                 else
%                     liftSurface(i) = wingtail(dihedral,semispan,chord,sweep,sections);
%                 end
                
                fixFactor = fixFactor/2;
                stayOut = -nHist;
            end
            
            % Sometimes straight unswept untapered wings cause issues
            % when merging, so ensure first partition sweep is non-zero
            if sweep(1) < 1
                
                sweep(1) = 1;
            end
            
            % Chord may get too low in "last resort" method, if so make all
            % chord values equal to half the aftbody length
            if chord(1) < 1
                
                chord(:) = parameters.AftL/2;
                parChord(:) = 0.5;
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
                        
                        semispan(1) = semispan(1) + 0.5*semispan(1);
                        parSemispan(1) = parSemispan(1) + 0.5*parSemispan(1);
                            
                        if Control
                            liftSurface(i) = wingtail(dihedral,semispan,chord,sweep,sections,control);
                        else
                            liftSurface(i) = wingtail(dihedral,semispan,chord,sweep,sections);
                        end
                end
                
            else
                
                % If too many attempts already, last resort reduce chord
                
                % Alter physical chord to be created, and particle chord to
                % be fed back to optimiser
                chord(1) = chord(1) - 0.1*chord(1);
                parChord(1) = parChord(1) - 0.1*parChord(1);
                
                % Sometimes straight unswept untapered wings cause issues
                % when merging, so ensure taper ratio is slightly < 1
                for j = 2:numel(chord)
                    if chord(j) >= chord(j-1)
                        chord(j) = 0.99*chord(j-1);
                        parChord(j) = 0.99*parChord(j-1);
                    end
                end
                
                % Re-create aerofoil
                if Control
                    liftSurface(i) = wingtail(dihedral,semispan,chord,sweep,sections,control);
                else
                    liftSurface(i) = wingtail(dihedral,semispan,chord,sweep,sections);
                end
                
            end
            
            % Re-try merge
            [liftSurface,aftbody,success,reason] = bodyfoil(aftbody,liftSurface,[xOffset,zOffset]);
            
            % Increase number of attempts and stayOut condition
            attempt = attempt + 1;
            stayOut = stayOut + 1;
            
            % Shift current reason down one in the history array to make 
            % room for new reason
            reasonHistory = circshift(reasonHistory,1);
            
            % Use to debug if configuration cannot be assembled
            if attempt > 90
               attempt
            end
            
        end
        
        %%
        liftSurface.Offset = [xOffset,zOffset];
        
        % No need to feedback if successful first time
        if attempt > 1
            % Feedback particle and physical positions to optimiser
            parPos(dihedralVar) = parDihedral;
            parPos(semispanVar) = parSemispan;
            parPos(chordVar) = parChord;
            parPos(sweepVar) = parSweep;
            parPos(xOffsetVar) = parxOffset;
            parPos(zOffsetVar) = parzOffset;

            configPos(dihedralVar) = dihedral;
            configPos(semispanVar) = semispan;
            configPos(chordVar) = chord;
            configPos(sweepVar) = sweep;
            configPos(xOffsetVar) = xOffset;
            configPos(zOffsetVar) = zOffset;
        end
        
        % Discretise aerofoils based on target length
        liftSurface = discwing(liftSurface);
        
        if Structure && i == 1
            
            wingbox.Location = [0.1,0.8]';
            wingbox.SparThick = 0.010;
            wingbox.SkinThick = 0.005;
            
            liftSurface = wingstructure(liftSurface,wingbox);
        end
        
    end
        
    
end

if Wing
    
    parameters.Chord = chord;
    parameters.Semispan = semispan;
    parameters.Sweep = sweep;
    parameters.Dihedral = dihedral;
    
    if Control
        parameters.Control = control;
    end
    
    % First liftSurface assumed to be wing
    wing = liftSurface(1);
    
    sumArea = sum(wing.Area);
    
    parameters.MAC = sum(wing.Area.*wing.MAC)/sumArea;
    parameters.Aref = sumArea*2;
    parameters.Wingspan = sum(wing.Span);
    
else
    
    liftSurface = '';
    
    parameters.Aref = aftbody.Area;
    parameters.MAC = aftbody.Length/2;parameters.BodyW = aftbody.Width;
    parameters.wingspan = [];
    
end

if Fore
    
        % Find radial aft body points for fore body interpolation
    % size-1 for panel based over
    radial = radialpoints(aftbody.Points);
    [dim_f,~] = size(radial);
    
    % Create nose and link together with aftbody by interpolating to
    % create forebody
    
    % If zero length or radius, nose must be pointed, therefore empty as
    % first forebody points will equal singular nose point
    if ~all([noseLength,noseRad])
        Nose = '';
        nosePoints.x = zeros(1,dim_f);
        nosePoints.y = zeros(1,dim_f);
        nosePoints.z = zeros(1,dim_f);
    else
        Nose = nose([noseRad,noseLength,zNoseOff],dim_f,noseForeLength);
        nosePoints = Nose.Points;
    end
    
    forebody = foreGen(foreBodyLength,nosePoints,radial,0.5);
    
end

%% Put assembly into single cell, delete any empty parts
assemblyProperties = {Nose,forebody,aftbody,liftSurface};
assemblyProperties(strcmp(assemblyProperties,'')) = [];