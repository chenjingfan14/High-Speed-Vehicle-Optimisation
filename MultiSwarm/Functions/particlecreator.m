function [assemblyProperties,Aref,MAC,flag] = particlecreator(configPos,partArrays,sections)
%% Assembles arbitrary configuration
% Will fuse together multiple aerofoil sections with body. Will also work
% for aerofoil alone or body alone configurations

%% Find where different parts exist within configuration
% Initialise existence of parts as false
[dim,~] = size(partArrays);

[wingCon,aftbodyCon,forebodyCon,noseCon] = deal(false(dim,1));

% Loop through configuration and set parts to true at that point in the 
% configuration cell when part name = "Wing" or "Nose" etc
for ii = 1:dim
    wingCon(ii) = strcmp(partArrays{ii,1},"Wing") | strcmp(partArrays{ii,1},"Tail");
    aftbodyCon(ii) = strcmp(partArrays{ii,1},"Aftbody");
    forebodyCon(ii) = strcmp(partArrays{ii,1},"Forebody");
    noseCon(ii) = strcmp(partArrays{ii,1},"Nose");
end

% Turn 1:dim logicals into indexes containing only numbers where said part
% type resides
array = 1:dim;
wingDim = array(wingCon);
aftbodyDim = array(aftbodyCon);
forebodyDim = array(forebodyCon);
noseDim = array(noseCon);

%% Create aerofoils
% If no aerofoils in configuration, loop will be skipped as dim = 0
% Count back from number of aerofoils for preallocation purposes
for ii=wingDim:-1:1
    
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
    
    % Create aerofoil
    aerofoil(ii) = wingtail(dihedral,semispan,chord,sweep,offset,sections);
end

%% Create body and assemble
% Check if body is part of configuration, if so build it, if not set body
% parts to empty. Currently aftbody defines body existence, and if it does
% nose and forebody must also exist

if any(aftbodyCon)
    
    % Again initialise correct indices along with their physical defition
    % for each body part
    % aftbodyStr = partArrays{aftbodyDim,2};
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
    xNosePlusFore = noseLength + foreBodyLength;
    
    % Create aft body portion of body
    aftbody = arbitraryfuse(aftbodyPos,xNosePlusFore);
    
    if ~any(wingCon)
        MAC = aftbody.Length/2;
        Aref = aftbody.Area;
        aerofoil = '';
    else
        % Join body and aerofoil together
        [aerofoil,aftbody,success] = bodyfoil(aftbody,aerofoil);
        
        % If joining is unsuccessful, configuration not feasible, stop
        % bodyfoil and return to particlecreator with flag
        if success
            flag = false;
        else
            flag = true;
            [assemblyProperties,Aref,MAC] = deal([]);
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
        Nose = nose(nosebodyPos,dim_f,xNosePlusFore);
        nosePoints = Nose.Points;
    end
    
    forebody = foreGen(forebodyPos,nosePoints,radial,0.5);
    
else % Set all body parts to empty
    aftbody = ''; forebody = ''; Nose = '';
end

%% New wing discretisation
if any(wingCon) % Aerofoil(s)
    
    % Discretise aerofoils based on target length
    aerofoil = discwing(aerofoil);
    % Aerofoil 1 assumed to be wing
    MAC = aerofoil(1).MAC;
    Aref = aerofoil(1).Area;
    
end


%% Put assembly into single cell, delete any empty parts
assemblyProperties = {Nose,forebody,aftbody,aerofoil};
assemblyProperties(strcmp(assemblyProperties,'')) = [];