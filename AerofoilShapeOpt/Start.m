%% Start Optimisation

% Initialise problem specific parameters, adds them to overall workspace
problemdefinition();

method = "PSO";
foils = getaerofoilsecdata;

tic

for i = 1:length(foils)
    
    options.Compare = foils{i};
    dim = ceil(length(foils{i})/2);
    options.UpperDisc = foils{i}(1:dim,1);
    options.LowerDisc = foils{i}(dim:end,1);
    
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
    
    [singleOpt,~] = size(GlobalBestPos);
    
    if aerofoilMethod == "BP3434"
        
        optFoil = BP3434(GlobalBestPos(singleOpt,:),varArray,wingPartitions,1,options.LowerDisc,options.UpperDisc);
        
    elseif aerofoilMethod == "Bezier"
        
        optFoil = Bezier3(GlobalBestPos(singleOpt,:),wingPartitions,options.BezierControlPoints,1,options.LowerDisc,options.UpperDisc);
        
    elseif aerofoilMethod == "BezierTC"
        
        [optFoil,cLine,tLine] = BezierTC(GlobalBestPos(singleOpt,:),wingPartitions,options.BezierControlPoints,1,options.LowerDisc,options.UpperDisc);
    end
    
    figure
    hold on
    axis equal
    plot(foils{i}(:,1),foils{i}(:,2),'k--');
    plot(optFoil{1}(:,1),optFoil{1}(:,2),'k');
    optStr = "Target";
    optStr(2) = ['Optimised ' aerofoilMethod];
    
    if exist('cLine','var')
        
        plot(cLine{1}(:,1),cLine{1}(:,2),'r');
        optStr(3) = "Camber";
    end
    if exist('tLine','var')
        
        plot(tLine{1}(:,1),tLine{1}(:,2),'b');
        optStr(4) = "Thickness";
    end
    legend(optStr)
    hold off
    
    iterations(i,:) = history(end,1);
    
    str = sprintf('AerofoilComparison%i', i);
    saveas(gcf,[resultPath, str]);
    
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