%% Start Particle Swarm Optimiser

% Initialise problem specific parameters, adds them to overall workspace
initialise();

tic

% Swarm size (must be divisible by 2 & 3 for global best and mutation
% subsets in MOPSO)
nPop = 120;
maxIt = 100; % Maximum number of iterations

w = 0.3; % Intertia coeff
c1 = 1.49; % Personal acceleration coeff
c2 = 1.49; % Social acceleration coeff

fi = maxIt; % Display Pareto Front evey fi iterations

%% Main PSO program
if nFun == 1 % Use Single-Objective Algorithm
    
    % If opt cost function value = max(f(x)) (rather than min) then can use
    % this to invert CF values for display purposes
    inv = false;
    % Max and min inertial values
    wmax = 0.8;
    wmin = 0.1;
    % Max stall values before simulation ends
    maxStall = 50;
    
    % Simplest Single-Objective PSO Algorithm
    % [GlobalBestFit,GlobalBestPos,history] = PSO(cond,costFun,varArray,varMin,varMax,nVar,nPop,maxIt,maxStall,w,wmax,wmin,c1,c2,nFun,inv,fi,foilData,n,flow,thetaBetaM,maxThetaBetaM,PrandtlMeyer,options);
    
    % Smarter Single-Objective PSO Algorithm
    [GlobalBestFit,GlobalBestPos,history] = PSONeighbourhood(cond,costFun,varArray,varMin,varMax,nVar,nPop,maxIt,maxStall,w,wmax,wmin,c1,c2,nFun,inv,fi,foilData,n,flow,thetaBetaM,maxThetaBetaM,PrandtlMeyer,options);

else % Use Multi-Objective Algorithm
    
    inv = false(1,nFun);
    maxPF = 100; % Maximum number of Pareto Front values
    mutProb = 1/nVar; % Probability of mutation
    
    [GlobalBestFit,GlobalBestPos,history] = MOPSO(cond,costFun,varArray,varMin,varMax,nVar,nPop,maxIt,maxPF,mutProb,w,c1,c2,nFun,inv,fi,foilData,n,flow,thetaBetaM,maxThetaBetaM,PrandtlMeyer,options);
end

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