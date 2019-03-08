function [Base] = baselinefun(variCons,flow,options,thetaBetaM,maxThetaBetaM,PrandtlMeyer,display)
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

wing = options.Wing;
aft = options.Aft;
fore = options.Fore;
nose = options.Nose;
control = options.Control;

Definition = {...
    "Variables",    "Values"};

if wing
    
    wingDefs = baseline_wing();
    Definition = [Definition; wingDefs];
end

if any([aft,fore,nose])
    
    baseline_body();
    Definition = [Definition; bodyDefs];
end

if control
    
    controlDefs = baseline_control();
    Definition = [Definition; controlDefs]; 
end

[foilData,~] = getaerofoilsecdata();

% If baseline is a direct definition (ie. actual physical values not
% needing to be transformed or conditioned) then find transformations in
% optimisation definition and invert variables so that they are of the same
% format as optimisation values
if isDirect
    
    [~,col] = size(variCons);
    
    % Grab titles
    for i=col:-1:1

        Headers(i) = variCons{1,i};
    end

    transCol = Headers == "Transformations";
    Definition = [Definition, variCons(:,transCol)];
end

[cond,varArray,baseVar] = translateOpt(Definition);

sectionPos = baseVar(any(varArray == ["Section","Bezier"],2));
sections = foilData(sectionPos);

baseVar = baseVar';

[baseProperties,allPoints,~,parameters] = particlecreator(baseVar,baseVar,varArray,sections,options);
% parameters.Aref = 33.213;

[~,Base.Results] = aeroprediction(baseProperties,allPoints,flow,parameters,thetaBetaM,maxThetaBetaM,PrandtlMeyer,options);

% Used to save direct inputs for postprocess
if isDirect
    baseVar = hardtransform(baseVar,cond,varArray,"inverse");
    
    % Inverse transform means that configuration is no longer direct
    % (necessary to state this for plotting)
    isDirect = false;
end

Base.Definition = cond;
Base.VarArray = varArray;
Base.Variables = baseVar;
Base.nVar = length(baseVar);
Base.Direct = isDirect;
Base.Sections = sections;

if wing
    
    Base.Bezier = options.Bezier;
    Base.nPartitions = length(parameters.Semispan);
end

if ~exist('display','var')
    
    display = false;
end

if display
    
    points = flowfinder(baseProperties);
    plotter(points,"title",'Baseline Configuration')
else
    save(fullfile([pwd '\Results'], 'Baseline'))
end
