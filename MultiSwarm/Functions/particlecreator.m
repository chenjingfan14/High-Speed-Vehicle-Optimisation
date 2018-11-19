function [assemblyProperties,parPos,configPos,parameters,success] = particlecreator(parPos,configPos,partArrays,sections,attempt)
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
    
    % Assign variables to their physical attribute
    dihedral = wingPos(wingStr == "Dihedral");
    chord = wingPos(wingStr == "Chord");
    sweep = wingPos(wingStr == "LESweep");
    semispan = wingPos(wingStr == "Semispan");
    offset = wingPos(any(wingStr == ["xOffset","zOffset"]',1));
    
    if attempt > 1
        
        chordID = partArrays{wingDim,2} == "Chord";
        ID = find(chordID,1,'first');
        
        % Alter physical chord to be created
        Cr = chord(1);
        Cr = Cr - 0.1*Cr;
        chord(1) = Cr;
        
        % Alter particle chord, feedback both to costcaller
        parChord = parPos(chordID);
        
        parCr = parChord(1);
        parCr = parCr - 0.1*parCr;
        parChord(1) = parCr;
        
        % Enforce taper <= 1
        for j = 2:numel(chord)
             if chord(j) > chord(j-1)
                 chord(j) = chord(j-1);
                 parChord(j) = parChord(j-1);
             end
        end
        
        parPos(chordID) = parChord;
        configPos(chordID) = chord;
        
    end
    
    % Create aerofoil
    if isempty(controlDim)
        liftSurface(i) = wingtail(dihedral,semispan,chord,sweep,offset,sections);
    else
        controlArray = partArrays{controlDim,3};
        control = configPos(controlArray);
        liftSurface(i) = wingtail(dihedral,semispan,chord,sweep,offset,sections,control);
    end
    
end

%% Create body and assemble
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
    
    if any(wingCon)
        % Join body and aerofoil together
        [liftSurface,aftbody,success] = bodyfoil(aftbody,liftSurface);
        
        % If joining is unsuccessful, configuration not feasible, stop
        % bodyfoil and return to particlecreator with flag
        if ~success
            [assemblyProperties,parameters] = deal([]);
            return
        end
    end
    
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
    
else % Set all body parts to empty
    aftbody = ''; forebody = ''; Nose = '';
end

%% New wing discretisation
if any(wingCon) % Aerofoil(s)
%     surfaceParameters = [0.5,0.8,0.7];
%     
%     liftSurface = controlsurface(liftSurface,surfaceParameters);
    
    % Discretise aerofoils based on target length
    liftSurface = discwing(liftSurface);
    % liftSurface 1 assumed to be wing
    
    wing = liftSurface(1);
    
    Aref = sum(wing.Area);
    MAC = sum(wing.Area.*wing.MAC)/Aref;
    
    wingspan = sum(wing.Span);
else
    Aref = aftbody.Area;
    MAC = aftbody.Length/2;parameters.BodyW = aftbody.Width;
    wingspan = [];
    liftSurface = '';
end

parameters.Aref = Aref;
parameters.Wingspan = wingspan;
parameters.MAC = MAC;
parameters.NoseL = noseForeLength;
parameters.BodyL = noseForeLength + aftBodyLength;
parameters.BodyW = aftbody.Width;
parameters.BodyH = aftbody.Height;

%% Put assembly into single cell, delete any empty parts
assemblyProperties = {Nose,forebody,aftbody,liftSurface};
assemblyProperties(strcmp(assemblyProperties,'')) = [];