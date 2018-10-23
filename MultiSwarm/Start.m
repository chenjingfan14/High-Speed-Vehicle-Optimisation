% Initialiser
clear all
close all
clc

options = simOptions();

% Swarm size (must be divisible by 3 for mutation subsets in MOPSO)
nPop = 300;
maxIt = 2000; % Maximum number of iterations

w = 0.3; % Intertia coeff
c1 = 1.49; % Personal acceleration coeff
c2 = 1.49; % Social acceleration coeff

% Number of decision variables (cost function values)
nFun = 4;

fi = maxIt; % Display Pareto Front evey fi iterations

% Cell containing variable names and their subsequent conditions to be 
% applied, so far these conditions are as follows:
% minimum value (mini) - if value is below those thresholds, set to zero
% floor - round value down to nearest integer
% if all - if any values = 0, set all to zero

% Number of wing partitions
n = 3;

variCons = {"Variables",    "Num Of",   "Conditions"    'Transformations';...
    "Dihedral",             "~",        "~",            '~';...
    "Chord",                n+1,        "< Previous",   '.*AftLength';...
    "LESweep",              n,          "~"             '~';...
    "Semispan",             n,          "Minimum 0.5",  '~';...            
    "Section",              n+1,        "Floor",        '~';...
    "xOffset",              "~",        "~",            '.*AftLength';...
    "zOffset",              "~",        "~"             '.*AftHeight/2';...
    "UpperLength",          "~",        "~",            '~';...
    "yUpperRad",            "~",        "~",            '~';...
    "yBotRatio",            "~",        "~",            '~';...
    "zUpperRad",            "~",        "~",            '~';...
    "SideLength",           "~",        "~",            '~';...
    "zLowerRad",            "~",        "~",            '~';...
    "AftLength",            "~",        "~",            '~';...
    "NoseRad",              "~",        "~",            '~';...
    "NoseLength",           "~",        "~",            '.*NoseRad';...
    "zNoseOffset",          "~",        "~",            '.*AftHeight/2';...
    "ForeLength",           "~",        "~",            '~'};
% Streamlines above cell to indices instead of names
% Total number of variables
[cond,varArray,nVar] = translate(variCons);

% Minimum and maximum variable values. Corresponds to variCons
%% Full Configuration (Wing & Body)
% varMin = [0, 0.5,0.1,0.1,0.1, 0,0,0, 2,0,0, 1,1,1,1, 0,-0.5,... % Wing
%     0,0.1,0.1, 0.1,0,0.1, 4, 0,0,-0.5, 1]; % Body
% varMax = [20, 1,1,1,1, 45,45,45, 5,2,2, 33,33,33,33, 0.5,0,...
%     1,1,1, 1,1,1, 10, 0.5,1,0, 5];

%% Wing Only
% varMin = [0, 0.5,0.1,0.1,0.1, 0,0,0, 2,0,0, 1,1,1,1, 0,-0.5,... % Wing
%     NaN,NaN,NaN, NaN,NaN,NaN, NaN, NaN,NaN,NaN, NaN]; % Body
% varMax = [20, 1,1,1,1, 45,45,45, 5,2,2, 33,33,33,33, 0.5,0,...
%     NaN,NaN,NaN, NaN,NaN,NaN, NaN, NaN,NaN,NaN, NaN];

%% Body Only
% varMin = [NaN, NaN,NaN,NaN,NaN, NaN,NaN,NaN, NaN,NaN,NaN, NaN,NaN,NaN,NaN, NaN,NaN,... % Wing
%     0,0.1,0.1, 0.1,0,0.1, 4, 0,0,-0.5, 1]; % Body
% varMax = [NaN, NaN,NaN,NaN,NaN, NaN,NaN,NaN, NaN,NaN,NaN, NaN,NaN,NaN,NaN, NaN,NaN,...
%     1,1,1, 1,1,1, 10, 0.5,1,0, 5];

%%
if exist('varMin','var')
    standard = isnan(varMin);
    [varMin,varMax] = standardvariables(standard,n,varMin,varMax);
else
    standard = true(1,nVar);
    [varMin,varMax] = standardvariables(standard,n);
end

% Load coordinates of 2D aerofoil sections into matrices within cell & 
% freestream flow parameters into structure
foilData = getaerofoilsecdata();

% Call flowparameters() to apply pre-defined angle attack/Mach numbers
% Call flowparameters(10,3) for example to run configurations at angle of
% attack 10deg and Mach 3. Multiple values can be input eg.
% flowparameters([0,2,4,6],[3,4]). See inside function for predetermined
% figures
flow = flowparameters();

% Load lookup tables for shock-expansion and Prandtl Meyer expansion
load('thetaBetaCurves.mat');
Mrange = [1:0.0001:10,10.1:0.1:100];
PrandtlMeyer = prandtlmeyerlookup(Mrange,flow);

costFun = @aeroprediction; % Cost function caller

% Main PSO program
if nFun == 1
    
    % If opt cost function value = max(f(x)) (rather than min) then can use
    % this to invert CF values for display purposes
    inv = false;
    % Max and min inertial values
    wmax = 0.8;
    wmin = 0.1;
    % Max stall values before simulation ends
    maxStall = 500;
    
    [GlobalBestFit,GlobalBestPos] = PSO(cond,costFun,varArray,varMin,varMax,nVar,nPop,maxIt,maxStall,w,wmax,wmin,c1,c2,nFun,inv,fi,foilData,flow,thetaBetaM,maxThetaBetaM,PrandtlMeyer,options);
else
    
    inv = false(1,nFun);
    inv = logical(inv);
    maxPF = nPop; % Maximum number of Pareto Front values
    mutProb = 1/nVar; % Probability of mutation

    [GlobalBestFit,GlobalBestPos] = MOPSO(cond,costFun,varArray,varMin,varMax,nVar,nPop,maxIt,maxPF,mutProb,w,c1,c2,nFun,inv,fi,foilData,flow,thetaBetaM,maxThetaBetaM,PrandtlMeyer,options);
end

% Use this function to create output plots of configurations
viewcaller(GlobalBestPos,cond,varArray,foilData,flow);