function [GlobalBestFit,GlobalBestPos,history] = GeneticAlgorithm(cond,costFun,varArray,varMin,varMax,nVar,nPop,maxIt,maxStall,Pc,Pm,Er,nFun,inv,fi,foilData,n,flow,thetaBetaM,maxThetaBetaM,PrandtlMeyer,options)

%%  Initialization

varMinMat = repmat(varMin',nPop,1);
varMaxMat = repmat(varMax',nPop,1);

varSize = [nPop nVar];

history = zeros(maxIt, nFun + 2);

%% Initialise

% Array to hold best cost value on each iteration
GlobalBestFit = zeros(maxIt+1,1);
GlobalBestPos = zeros(maxIt+1,nVar);

% Create population array
popGene = unifrnd(varMinMat,varMaxMat,varSize);

[popFitness,popGene] = costcaller(costFun,nPop,nFun,popGene,cond,varArray,n,foilData,flow,thetaBetaM,maxThetaBetaM,PrandtlMeyer,options);

% for i = 1 : M
%     population.Chromosomes(i).fitness = Problem.obj( population.Chromosomes(i).Gene(:) );
% end

[GlobalBestFit(1),bestID] = min(popFitness);
GlobalBestPos(1,:) = popGene(bestID,:);

GlobalBestFitDisp = GlobalBestFit;
GlobalBestFitDisp(:,inv) = inv(:,inv)./GlobalBestFitDisp(:,inv);

history(1,:) = [1, GlobalBestFitDisp(1,:), 0];

fitnessBar = mean(popFitness,1);
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
        
        % Mutation
        c1Gene = mutation_continious(c1Gene, Pm, varMin, varMax);
        c2Gene = mutation_continious(c2Gene, Pm, varMin, varMax);
        
        newPopGene([k-1 k],:) = [c1Gene; c2Gene];
    end
        
    [newPopFitness,newPopGene] = costcaller(costFun,nPop,nFun,newPopGene,cond,varArray,n,foilData,flow,thetaBetaM,maxThetaBetaM,PrandtlMeyer,options);
    
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
    
    if GlobalBestFit(it) < GlobalBestFit(it-1)

        stall = 0; 
    else
        stall = stall+1;
    end
    
    fitnessBar = mean(popFitness,1);
    fprintf('Iteration %i: Global Best: %3.4f Mean: %3.4f Stall: %i \n', it-1, GlobalBestFitDisp(it), fitnessBar, stall);
    
    history(it,:) = [it-1, GlobalBestFitDisp(it), stall];
    
    if mod(it-1,fi)==0 || stall == maxStall
        
        figure(fcount)
        clf
        title(['Iteration: ' num2str(it-1)])
        set(gcf, 'Position', [0, 0, 1920, 1200])
        hold on
        xlabel('Iterations')
        ylabel('f(x)')
        plot(1:it,GlobalBestFitDisp(1:it));
        hold off
        % Pause to display graph while simulation is running
        pause(0.00001)
        fcount=fcount+1;
        
        if stall == maxStall
            
            GlobalBestFit(it+1:end) = [];
            GlobalBestPos(it+1:end,:) = [];
            history(it+1:end,:) = [];
            break
        end
    end
end

end