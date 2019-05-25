%% Start Particle Swarm Optimiser

tic

% Swarm size (must be divisible by 2 & 3 for global best and mutation
% subsets in MOPSO)
nPop = 120;
maxIt = 100; % Maximum number of iterations

w = 0.3; % Intertia coeff
c1 = 1.49; % Personal acceleration coeff
c2 = 1.49; % Social acceleration coeff
mutProb = 1/nVar; % Probability of mutation

fi = maxIt; % Display Pareto Front evey fi iterations

%% Main PSO program
if nFun == 1 % Use Single-Objective Algorithm
    
    % Max and min inertial values
    wmax = 0.8;
    wmin = 0.1;
    % Max stall values before simulation ends
    maxStall = 50;
    
    % Simplest Single-Objective PSO Algorithm
    % [GlobalBestFit,GlobalBestPos,history] = PSO(cond,costFun,varArray,varMin,varMax,nVar,nPop,maxIt,maxStall,mutProb,w,wmax,wmin,c1,c2,nFun,inv,neg,fi,foilData,n,flow,thetaBetaM,maxThetaBetaM,PrandtlMeyer,options);
    
    % Smarter Single-Objective PSO Algorithm
    [GlobalBestFit,GlobalBestPos,history] = PSONeighbourhood(cond,costFun,varArray,varMin,varMax,nVar,nPop,maxIt,maxStall,mutProb,w,wmax,wmin,c1,c2,nFun,fi,options);
    
else % Use Multi-Objective Algorithm

    maxPF = 100; % Maximum number of Pareto Front values
    
    [GlobalBestFit,GlobalBestPos,history] = MOPSO(cond,costFun,varArray,varMin,varMax,nVar,nPop,maxIt,maxPF,mutProb,w,c1,c2,nFun,fi,options);
end

time = toc;