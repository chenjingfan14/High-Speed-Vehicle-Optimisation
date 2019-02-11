%% Start Genetic Algorithm

% Initialise problem specific parameters, adds them to overall workspace
initialise();

tic

addpath(genpath('GeneticAlgorithm'))

% Swarm size (must be divisible by 2 & 3 for global best and mutation
% subsets in MOPSO)
nPop = 60;
maxIt = 50; % Maximum number of iterations

w = 0.3; % Intertia coeff
c1 = 1.49; % Personal acceleration coeff
c2 = 1.49; % Social acceleration coeff

% Number of decision variables (cost function values)
nFun = 1;

fi = maxIt; % Display Pareto Front evey fi iterations

%% Main GA program (Only single-objective currently available)

% If opt cost function value = max(f(x)) (rather than min) then can use
% this to invert CF values for display purposes
inv = false;
% Max stall values before simulation ends
maxStall = 50;

%% Controling paramters of the GA algortihm
Pc = 0.95;         % Probablility of crossover
Pm = 1/nVar;        % Probability of mutation 
Er = 0.2;          % Elitism ratio 

[GlobalBestFit,GlobalBestPos,history] = GeneticAlgorithm(cond,costFun,varArray,varMin,varMax,nVar,nPop,maxIt,maxStall,Pc,Pm,Er,nFun,inv,fi,foilData,n,flow,thetaBetaM,maxThetaBetaM,PrandtlMeyer,options);

% If running on cluster, close down parallel loop
if cluster
    delete(gcp('nocreate'));
end

time = toc;

% configInputs variable required for postprocess function
configInputs = GlobalBestPos;

% Save global best history/pareto front figure and workspace in current directory
if nFun == 1
    
    saveas(gcf,'GlobalBestHistory')
else
    saveas(gcf,'ParetoFront')
end

save('OptimisationResults')

% Use this function to create output plots/results for arbitrary configs, 
% will not work for cluster simulations
if ~cluster
    postprocess
end