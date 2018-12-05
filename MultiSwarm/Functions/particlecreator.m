function [assemblyProperties,parPos,parameters] = particlecreator(parPos,configPos,partArrays,sections)
%% Assembles arbitrary configuration
% Will fuse together multiple aerofoil sections with body. Will also work
% for aerofoil alone or body alone configurations

%% Find where different parts exist within configuration
% Initialise existence of parts as false
[dim,~] = size(partArrays);

[wingCon,aftbodyCon,forebodyCon,noseCon,controlCon] = deal(false(dim,1));

% Loop through configuration and set parts to true at that point in the
% configuration cell when part name = "Wing" or "Nose" etc
for i = 1:dim
    wingCon(i) = strcmp(partArrays{i,1},"Wing") | strcmp(partArrays{i,1},"Tail");
    aftbodyCon(i) = strcmp(partArrays{i,1},"Aftbody");
    forebodyCon(i) = strcmp(partArrays{i,1},"Forebody");
    noseCon(i) = strcmp(partArrays{i,1},"Nose");
    controlCon(i) = strcmp(partArrays{i,1},"Control");
end

% Turn 1:dim logicals into indexes containing only numbers where said part
% type resides
array = 1:dim;
wingDim = array(wingCon);
aftbodyDim = array(aftbodyCon);
forebodyDim = array(forebodyCon);
noseDim = array(noseCon);
controlDim = array(controlCon);

%% Create body
% Check if body is part of configuration, if so build it, if not set body
% parts to empty. Currently aftbody defines body existence, and if it does
% nose and forebody must also exist

if any(aftbodyCon)
    
    % Again initialise correct indices along with their physical defition
    % for each body part
    aftbodyStr = partArrays{aftbodyDim,2};
    forebodyStr = partArrays{forebodyDim,2};
    noseStr = partArrays{noseDim,2};
    
    aftbodyArray = partArrays{aftbodyDim,3};
    forebodyArray = partArrays{forebodyDim,3};
    noseArray = partArrays{noseDim,3};
    
    aftbodyPos = configPos(aftbodyArray);
    forebodyPos = configPos(forebodyArray);
    nosebodyPos = configPos(noseArray);
    
    noseLength = nosebodyPos(noseStr == "NoseLength");
    noseRad = nosebodyPos(noseStr == "NoseRad");
    foreBodyLength = forebodyPos(forebodyStr == "ForeLength");
    noseForeLength = noseLength + foreBodyLength;
    
    aftBodyLength = aftbodyPos(aftbodyStr == "AftLength");
    % Create aft body portion of body
    aftbody = arbitraryfuse(aftbodyPos,noseForeLength);
    
    parameters.NoseL = noseForeLength;
    parameters.BodyL = noseForeLength + aftBodyLength;
    parameters.BodyW = aftbody.Width;
    parameters.BodyH = aftbody.Height;
    
else % Set all body parts to empty
    aftbody = ''; forebody = ''; Nose = '';
end

%% Create aerofoils
% If no aerofoils in configuration, loop will be skipped as dim = 0
% Count back from number of aerofoils for preallocation purposes
for i=wingDim:-1:1
    
    % Initialise what indices of configPos refer to wing, along with what
    % each position physcially is (eg. "Chord", "Semispan" etc)
    wingStr = partArrays{wingDim,2};
    wingArray = partArrays{wingDim,3};
    
    % Use above indices to take only wing related parameters for configPos
    wingPos = configPos(wingArray);
    wingPar = parPos(wingArray);
    
    dihedralID = wingStr == "Dihedral";
    chordID = wingStr == "Chord";
    sweepID = wingStr == "LESweep";
    semispanID = wingStr == "Semispan";
    xOffsetID = wingStr == "xOffset";
    zOffsetID = wingStr == "zOffset";
    
    % Assign variables to their physical attribute
    dihedral = wingPos(dihedralID);
    chord = wingPos(chordID);
    sweep = wingPos(sweepID);
    semispan = wingPos(semispanID);
    xOffset = wingPos(xOffsetID);
    zOffset = wingPos(zOffsetID);
    
    parDihedral = wingPar(dihedralID);
    parChord = wingPar(chordID);
    parSweep = wingPar(sweepID);
    parSemispan = wingPar(semispanID);
    parxOffset = wingPar(xOffsetID);
    parzOffset = wingPar(zOffsetID);
    
    % If first semispan = 0, give starting value so that wing can be
    % created
    if semispan == 0
        semispan(1) = 0.1;
        parSemispan(1) = 0.1;
    end
    
    % Create aerofoil
    if isempty(controlDim)
        liftSurface(i) = wingtail(dihedral,semispan,chord,sweep,sections);
    else
        controlArray = partArrays{controlDim,3};
        control = configPos(controlArray);
        liftSurface(i) = wingtail(dihedral,semispan,chord,sweep,sections,control);
    end
    
    if any(aftbodyCon)
        
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
            
            % Sometimes straight unswept untapered wings cause issues
            % when merging, so ensure first partition sweep is non-zero
            if sweep(1) < 1
                sweep(1) = 1;
            end
            
            % Chord may get too low in "last resort" method, if so make all
            % chord values equal to half the entire body length
            if chord(1) < 1
                chord(:) = parameters.BodyL/2;
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
                        
                    case 5
                        
                        % If inboard semispan = 0, it will not increase
                        % below, thus add small number to make it non-zero
                        if semispan(1) == 0
                            semispan(1) = 0.1;
                            parSemispan(1) = 0.1;
                        else % Increase semispan
                            semispan(1) = semispan(1) + 0.1*semispan(1);
                            parSemispan(1) = parSemispan(1) + 0.1*parSemispan(1);
                        end
                            
                        if isempty(controlDim)
                            liftSurface(i) = wingtail(dihedral,semispan,chord,sweep,sections);
                        else
                            liftSurface(i) = wingtail(dihedral,semispan,chord,sweep,sections,control);
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
                if isempty(controlDim)
                    liftSurface(i) = wingtail(dihedral,semispan,chord,sweep,sections);
                else
                    liftSurface(i) = wingtail(dihedral,semispan,chord,sweep,sections,control);
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
        
        % Feedback particle and physical positions to optimiser
        parPos(wingArray(dihedralID)) = parDihedral;
        parPos(wingArray(semispanID)) = parSemispan;
        parPos(wingArray(chordID)) = parChord;
        parPos(wingArray(sweepID)) = parSweep;
        parPos(wingArray(xOffsetID)) = parxOffset;
        parPos(wingArray(zOffsetID)) = parzOffset;

        configPos(wingArray(dihedralID)) = dihedral;
        configPos(wingArray(semispanID)) = semispan;
        configPos(wingArray(chordID)) = chord;
        configPos(wingArray(sweepID)) = sweep;
        configPos(wingArray(xOffsetID)) = xOffset;
        configPos(wingArray(zOffsetID)) = zOffset;
        
        % Discretise aerofoils based on target length
        liftSurface = discwing(liftSurface);
        
    end
        
    
end

if any(wingCon)
    
    parameters.Chord = chord;
    parameters.Semispan = semispan;
    parameters.Sweep = sweep;
    parameters.Dihedral = dihedral;
    
    if ~isempty(controlDim)
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

if any(forebodyCon)
    
        % Find radial aft body points for fore body interpolation
    % size-1 for panel based over
    radial = radialpoints(aftbody.Points);
    [~,dim_f] = size(radial.y);
    
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
        Nose = nose(nosebodyPos,dim_f,noseForeLength);
        nosePoints = Nose.Points;
    end
    
    forebody = foreGen(forebodyPos,nosePoints,radial,0.5);
    
end

%% Put assembly into single cell, delete any empty parts
assemblyProperties = {Nose,forebody,aftbody,liftSurface};
assemblyProperties(strcmp(assemblyProperties,'')) = [];