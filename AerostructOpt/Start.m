%% Start Optimisation

% Initialise problem specific parameters, adds them to overall workspace
problemdefinition();

tic

% Population size (must be divisible by 2 & 3 for global best and mutation
% subsets in MOPSO)
nPop = 120;
maxIt = 100; % Maximum number of iterations
Pm = 1/nVar; % Probability of mutation

fi = maxIt; % Display Pareto Front evey fi iterations
        
% Begin Optimisation
switch method
    
    case "PSO"
        
        w = 0.3; % Intertia coeff
        c1 = 1.49; % Personal acceleration coeff
        c2 = 1.49; % Social acceleration coeff
        
        if nFun == 1 % Use Single-Objective Algorithm
            
            % Max and min inertial values
            wmax = 0.8;
            wmin = 0.1;
            % Max stall values before simulation ends
            maxStall = 50;
            
            % Simplest Single-Objective PSO Algorithm
            % [GlobalBestFit,GlobalBestPos,history] = PSO(cond,costFun,varArray,varMin,varMax,nVar,nPop,maxIt,maxStall,mutProb,w,wmax,wmin,c1,c2,nFun,inv,neg,fi,foilData,n,flow,thetaBetaM,maxThetaBetaM,PrandtlMeyer,options);
            
            % Smarter Single-Objective PSO Algorithm
            [GlobalBestFit,GlobalBestPos,history] = PSONeighbourhood(cond,costFun,varArray,varMin,varMax,nVar,nPop,maxIt,maxStall,Pm,w,wmax,wmin,c1,c2,nFun,fi,options);
            
        else % Use Multi-Objective Algorithm
            
            maxPF = 100; % Maximum number of Pareto Front values
            
            [GlobalBestFit,GlobalBestPos,history] = MOPSO(cond,costFun,varArray,varMin,varMax,nVar,nPop,maxIt,maxPF,Pm,w,c1,c2,nFun,fi,options);
        end
        
    case "GA"
        
        % Only single-objective currently available
        
        % Max stall values before simulation ends
        maxStall = 50;
        
        % Controling paramters of the GA algortihm
        Pc = 0.95;      % Probablility of crossover
        Er = 0.2;       % Elitism ratio
        
        [GlobalBestFit,GlobalBestPos,history] = GeneticAlgorithm(cond,costFun,varArray,varMin,varMax,nVar,nPop,maxIt,maxStall,Pc,Pm,Er,nFun,fi,options);
end

time = toc;

% If running on cluster, close down parallel loop
if cluster
    
    delete(gcp('nocreate'));
end

% configInputs variable required for postprocess function
configInputs = GlobalBestPos;

% Save global best history/pareto front figure and workspace in current directory
if nFun == 1
    
    saveas(gcf,'GlobalBestHistory')
%     saveas(gcf,[resultPath '\GlobalBestHistory'])
else
    saveas(gcf,'ParetoFront')
%     saveas(gcf,[resultPath '\ParetoFront'])
end

save('OptimisationResults')
% save(fullfile(resultPath, 'OptimisationResults'))

% Use this function to create output plots/results for arbitrary configs,
% will not work for cluster simulations
if ~cluster
    
    postprocess
end