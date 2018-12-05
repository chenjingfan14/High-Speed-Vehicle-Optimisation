function [base] = baselinefun(flow,options,thetaBetaM,maxThetaBetaM,PrandtlMeyer)
%% Baseline configuration to be improved upon (currently X-34)

% X-34 build uses pre-defined aerofoil sections so set Bezier to false
options.Bezier = false;
% Baseline not yet created for cost function so set to false
options.baseline = false;

% Assumed that baseline configs will input direct values and not need to be
% transformed. Thus will have to be reverse transformed to use - 
% for example body baseline variables - as standard variables (those to be 
% held cosntant)
isDirect = true;

control = options.control;

% Number of wing partitions
n = 2;

% No conditions required. However if variables are being input directly
% (ie. not being transformed same as initialisation variCons) then inverse
% transforms must be done here
variCons = {"Variables",    "Num Of",  "Conditions",  'Transformations';...
    "Dihedral",             "~",       "~",           '~';...
    "Chord",                n+1,       "~",           './AftLength';...
    "LESweep",              n,         "~"            '~';...
    "Semispan",             n,         "~",           '~';...            
    "SectionDefinition",    "~",       "~",           '~';...
    "xOffset",              "~",       "~",           './AftLength';...
    "zOffset",              "~",       "~"            './AftHalfHeight';...
    "UpperLength",          "~",       "~",           '~';...
    "yUpperRad",            "~",       "~",           '~';...
    "yBotRatio",            "~",       "~",           '~';...
    "zUpperRad",            "~",       "~",           '~';...
    "SideLength",           "~",       "~",           '~';...
    "zLowerRad",            "~",       "~",           '~';...
    "AftLength",            "~",       "~",           '~';...
    "NoseRad",              "~",       "~",           '~';...
    "NoseLength",           "~",       "~",           './NoseRad';...
    "zNoseOffset",          "~",       "~",           './AftHalfHeight';...
    "ForeLength",           "~",       "~",           '~'};

for i=size(variCons,1):-1:1
    secDef(i) = variCons{i,1} == "SectionDefinition";
end

baseWing = [4, 13.85,4.387,1.5241, 80,45, 0.744+0.88,2.62686, 2,1,1, -3.5,-0.74];
baseBody = [0,0.88,0.1, 0.88,0.759,0.05, 11.8956, 0.155,0.1115,-0.6, 4.4238];

wingEnd = length(baseWing);
bodyBegin = wingEnd + 1;

baseVar = [baseWing, baseBody];

if control
   
    controlVar = [0.4,0.7,0.7];
    baseVar = [baseVar, controlVar];
    
    variCons(end+1,:) = {"ControlSpan", 2, "> Previous", '~'};
    variCons(end+1,:) = {"ControlChord", "~", "~", '~'};
    
end

SectionDefinition = {"Section", n+1, "~", '~'};
variCons(secDef,:) = SectionDefinition;

[foilData,~] = getaerofoilsecdata();

[cond,varArray,~] = translate(variCons);
[partArrays,sectionArray] = partIndexing(cond,varArray);

sectionPos = baseVar(:,sectionArray);
sections = foilData(sectionPos);

[baseProperties,~,parameters] = particlecreator(baseVar,baseVar,partArrays,sections);
% parameters.Aref = 33.213;

[~,base] = aeroprediction(baseProperties,flow,parameters,thetaBetaM,maxThetaBetaM,PrandtlMeyer,options);

% Used to save direct inputs for postprocess
if isDirect
    configInputs = baseVar;
end

% Reverse transform into ratios etc that will be used in simulation
baseVar = reversetransform(configInputs,cond,varArray);

base.Wing = baseVar(1:wingEnd);
base.Body = baseVar(bodyBegin:end);

% Save basline so that it can be loaded into postprocess
save('Baseline')