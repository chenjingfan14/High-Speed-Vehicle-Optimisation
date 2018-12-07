function [Base] = baselinefun(flow,options,thetaBetaM,maxThetaBetaM,PrandtlMeyer)
%% Baseline configuration to be improved upon (currently X-34)

% X-34 build uses pre-defined aerofoil sections so set Bezier to false
options.Bezier = false;
% Baseline not yet created for cost function so set to false
options.Baseline = false;

% Assumed that baseline configs will input direct values and not need to be
% transformed. Thus will have to be reverse transformed to use - 
% for example body baseline variables - as standard variables (those to be 
% held cosntant)
isDirect = true;

% No conditions required. However if variables are being input directly
% (ie. not being transformed same as initialisation variCons) then inverse
% transforms must be done here
Definition = {...
    "Variables",    "Values"};
    
wingDefs = {...
    "Dihedral",     4;
    "Chord",       [13.85, 4.387, 1.5241];
    "LESweep",     [80, 45];
    "Semispan",    [0.744+0.88, 2.62686];           
    "Section",     ["FSPLStrake", "FSPLND", "FSPLND"];
    "xOffset",     -3.5;
    "zOffset",     -0.74};

aftDefs = {...
    "UpperLength",  0;
    "yUpperRad",    0.88;
    "yBotRatio",    0.1;
    "zUpperRad",    0.88;
    "SideLength",   0.759;
    "zLowerRad",    0.05;
    "AftLength",    11.8956};

foreDefs = {...
    "ForeLength",   4.4238};

noseDefs = {...
    "NoseRad",      0.155;
    "NoseLength",   0.1115;
    "zNoseOffset", -0.6};

controlDefs = {...
    "ControlSpan", [0.4,0.7];
    "ControlChord", 0.7};

Wing = options.Wing;
Aft = options.Aft;
Fore = options.Fore;
Nose = options.Nose;
Control = options.Control;

if Wing
    Definition = [Definition; wingDefs];
end

if Aft
    Definition = [Definition; aftDefs];
end

if Fore
    Definition = [Definition; foreDefs];
end

if Nose
    Definition = [Definition; noseDefs];
end

if Control
    Definition = [Definition; controlDefs]; 
end


[foilData,~] = getaerofoilsecdata();

[cond,varArray,baseVar] = translateOpt(Definition);

sectionPos = baseVar(:,any(varArray == ["Section","Bezier"],2));
sections = foilData(sectionPos,2);

% [wingDim,aftbodyDim,forebodyDim,noseDim,controlDim] = findparts(partArrays);

[baseProperties,~,parameters] = particlecreator(baseVar,baseVar,varArray,sections,options);
% parameters.Aref = 33.213;

[~,Base.Results] = aeroprediction(baseProperties,flow,parameters,thetaBetaM,maxThetaBetaM,PrandtlMeyer,options);

% Used to save direct inputs for postprocess
if isDirect
    configInputs = baseVar;
end

% Reverse transform into ratios etc that will be used in simulation
baseVar = hardtransform(configInputs,cond,varArray);

Base.VarArray = varArray;
Base.Variables = baseVar;
Base.nVar = length(baseVar);
Base.Bezier = options.Bezier;
Base.nPartitions = n;

% Save basline so that it can be loaded into postprocess
save('Baseline')