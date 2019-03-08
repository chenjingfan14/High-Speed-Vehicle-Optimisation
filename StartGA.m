%% Start Genetic Algorithm

% Initialise problem specific parameters, adds them to overall workspace
initialise();

tic

%% Main GA program (Only single-objective currently available)

% nPop must be even
nPop = 120;
maxIt = 100; % Maximum number of iterations

fi = maxIt; % Display Pareto Front evey fi iterations

% If opt cost function value = max(f(x)) (rather than min) then can use
% this to invert CF values for display purposes
inv = false;
% Max stall values before simulation ends
maxStall = 50;

% Controling paramters of the GA algortihm
Pc = 0.95;      % Probablility of crossover
Pm = 1/nVar;    % Probability of mutation 
Er = 0.2;       % Elitism ratio 

[GlobalBestFit,GlobalBestPos,history] = GeneticAlgorithm(cond,costFun,varArray,varMin,varMax,nVar,nPop,maxIt,maxStall,Pc,Pm,Er,nFun,inv,fi,foilData,n,flow,thetaBetaM,maxThetaBetaM,PrandtlMeyer,options);

%% Begin post-processing
% If running on cluster, close down parallel loop
if cluster
    delete(gcp('nocreate'));
end

time = toc;

% configInputs variable required for postprocess function
configInputs = GlobalBestPos;

% Save global best history/pareto front figure and workspace in current directory
if nFun == 1
    
    saveas(gcf,[resultPath '\GlobalBestHistory'])
else
    saveas(gcf,[resultPath '\ParetoFront'])
end

save(fullfile(resultPath, 'OptimisationResults'))

% Use this function to create output plots/results for arbitrary configs, 
% will not work for cluster simulations
if ~cluster
    postprocess
end