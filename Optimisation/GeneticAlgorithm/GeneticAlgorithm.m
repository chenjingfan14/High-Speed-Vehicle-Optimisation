function [GlobalBestFit,GlobalBestPos,history] = GeneticAlgorithm(cond,costFun,varArray,varMin,varMax,nVar,nPop,maxIt,maxStall,Pc,Pm,Er,nFun,fi,options)

%%  Initialization

varMinMat = repmat(varMin',nPop,1);
varMaxMat = repmat(varMax',nPop,1);

varSize = [nPop nVar];

history = zeros(maxIt, nFun + 2);

%% Initialise

tolerance = options.Tolerance;

% Array to hold best cost value on each iteration
GlobalBestFit = zeros(maxIt+1,1);
GlobalBestPos = zeros(maxIt+1,nVar);

% Create population array
popGene = unifrnd(varMinMat,varMaxMat,varSize);

[popFitness,popGene] = costcaller(costFun,nPop,nFun,popGene,cond,varArray,options);

[GlobalBestFit(1),bestID] = min(popFitness);
GlobalBestPos(1,:) = popGene(bestID,:);

inv = options.Inv;
neg = options.Neg;

GlobalBestFitDisp = GlobalBestFit;
GlobalBestFitDisp(:,inv) = inv(:,inv)./GlobalBestFitDisp(:,inv);
GlobalBestFitDisp(:,neg) = -GlobalBestFitDisp(:,neg);

history(1,:) = [1, GlobalBestFitDisp(1,:), 0];

fitnessBar = mean(popFitness(popFitness < inf),1);
fprintf('Iteration %i: Global Best: %3.4f Mean: %3.4f Stall: %i \n', 0, GlobalBestFitDisp(1,:), fitnessBar, 0);

stall = 0;
fcount = 1;

%% Main loop
for it = 2 : maxIt + 1
    
    for k = nPop: -2: 1
        % Selection
        [p1Gene,p2Gene] = selection(nPop,popFitness,popGene);
        
        % Crossover
        [c1Gene,c2Gene] = crossover_continious(p1Gene, p2Gene, Pc, nVar, varMin, varMax);
        
        newPopGene([k-1 k],:) = [c1Gene; c2Gene];
    end
    
    % Uniform mutation
    newPopGene = unimutation(newPopGene,nPop,nVar,varMinMat,varMaxMat,Pm);
        
    [newPopFitness,newPopGene] = costcaller(costFun,nPop,nFun,newPopGene,cond,varArray,options);
    
    [newPopFitness,sortID] = sort(newPopFitness,'descend');
    newPopGene = newPopGene(sortID,:);
    
    % Elitism
    [newPopFitness,newPopGene] = elitism(nPop, popFitness, popGene, newPopFitness, newPopGene, Er);
    
    popGene = newPopGene;
    popFitness = newPopFitness;
    
    [GlobalBestFit(it) , bestID ] = min(popFitness);
    GlobalBestPos(it,:) = popGene(bestID,:);
    
    %% Display
    GlobalBestFitDisp = GlobalBestFit;
    GlobalBestFitDisp(:,inv) = inv(:,inv)./GlobalBestFitDisp(:,inv);
    GlobalBestFitDisp(:,neg) = -GlobalBestFitDisp(:,neg);
    
    if GlobalBestFit(it) < GlobalBestFit(it-1)

        stall = 0; 
    else
        stall = stall+1;
    end
    
    fitnessBar = mean(popFitness(popFitness < inf),1);
    fprintf('Iteration %i: Global Best: %3.4f Mean: %3.4f Stall: %i \n', it-1, GlobalBestFitDisp(it), fitnessBar, stall);
    
    history(it,:) = [it-1, GlobalBestFitDisp(it), stall];
    
    convergence = stall == maxStall || GlobalBestFit(it) < tolerance;
    
    if mod(it-1,fi) == 0 || convergence
        
        if options.Baseline
            
            baselineCost = options.Base.Results.Cost;
        else
            baselineCost = zeros(1,nFun);
        end
        
        figure(fcount)
        clf
        title(['Cost Function Time History (Iteration: ' num2str(it-1) ')'])
        set(gcf, 'Position', [0, 0, 1920, 1200])
        hold on
        xlabel('Iterations')
        ylabel('f(x)')
        plot(it-1 ,baselineCost,'kx');
        plot(0:it-1 ,GlobalBestFitDisp(1:it),'k');
        legend('Baseline', 'Optimal Design')
        hold off
        % Pause to display graph while simulation is running
        pause(0.00001)
        fcount = fcount+1;
        
        if convergence
            
            GlobalBestFit(it+1:end) = [];
            GlobalBestPos(it+1:end,:) = [];
            history(it+1:end,:) = [];
            break
        end
    end
end

end