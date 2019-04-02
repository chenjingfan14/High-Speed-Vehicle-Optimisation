function [nonDomFit,nonDomPos,history] = MOPSO(cond,costFun,varArray,varMin,varMax,nVar,nPop,maxIt,maxPF,mutProb,w,c1,c2,nFun,inv,fi,options)
%% Multi Objective Particle Swarm Optimiser
% Main program, initialises swarm based on minimum/maximum design variables
% Uses cost function values to update swarm positions throughout process

% Particle numbers
nPopArray = (1:nPop)';

% Creating 3 sets within population, defines mutation method on set
% However all particles still communicate with each other regardless of
% sets
setSize = nPop/3;
set1 = 1:setSize; % No mutation
set2 = set1 + setSize; % Uniform mutation
set3 = set2 + setSize; % Non-uniform mutation

% Replicating min/max bounds into matrix with number of rows = population
% size
varMinMat = repmat(varMin',nPop,1);
varMaxMat = repmat(varMax',nPop,1);

varSize = [nPop nVar];

history = zeros(maxIt, nFun + 2);

%% Initialise Swarm

% Create population, within set bounds matrix created above,
% set initial velocity to zero
parPos = unifrnd(varMinMat,varMaxMat,varSize);
parVel = zeros(varSize);

%% Calculate initial fitness functions
[parFit,parPos] = costcaller(costFun,nPop,nFun,parPos,cond,varArray,options);

%%
% Initialise Pareto Front matrices
j = 1;
parsDominated = zeros([maxPF 1]);
nonDomFit = zeros([maxPF nFun]);
nonDomPos = zeros([maxPF nVar]);

% Checks if particle is dominated by any other, non-dominated particles are
% added to the Pareto Front and ranked based on how many other particles 
% they dominate
for i=1:nPop
    contender = parFit(i,:);
    % If particle dominated by any other
    con1 = all(parFit <= contender,2);
    con2 = any(parFit(con1,:) < contender);
    
    % If no to above, particle is non-dominated
    if ~any(con2)
        % How many particles does non-dominated particle dominate
        con1 = all(contender <= parFit,2);
        con2 = any(contender < parFit(con1,:),2);
        parsDominated(j) = sum(con2);
        nonDomFit(j,:) = parFit(i,:);
        nonDomPos(j,:) = parPos(i,:);
        j = j + 1;
    end
end

% Remove empty slots
parsDominated(j:end) = [];
nonDomFit(j:end,:) = [];
nonDomPos(j:end,:) = [];

% Max and min non-dominated particle fitness functions 
maxf = max(nonDomFit,[],1);
minf = min(nonDomFit,[],1);

% Normalise for crowding distance calculation
nonDomFitNorm = (nonDomFit - minf)./(maxf - minf);

% Find crowding distance of non-dominated particles
nonDomCrowd = crowder(nonDomFitNorm,maxf,minf,nFun);

% If non-dominated particles > max Pareto Front, discard those which
% are most crowded

[nNonDomParticles,~] = size(nonDomCrowd);

if nNonDomParticles > maxPF
    dim = (1:maxPF)';
    num = (1:size(nonDomCrowd,1))';
    csort = sortrows([num,nonDomCrowd],2,'descend');
    
    % Only need this if gbest selector needs it
    % nonDomCrowd = csort(dim,2);
    
    nonDomFit = nonDomFit(csort(dim,1),:);
    nonDomPos = nonDomPos(csort(dim,1),:);
    parsDominated = parsDominated(csort(dim,1));
    
    nNonDomParticles = maxPF;
    
end

% Initial Pareto Front = initial non-dominated particles
initPF = nonDomFit;
initPFDisp = initPF;
initPFDisp(:,inv) = inv(:,inv)./initPFDisp(:,inv);

% Initial particle best position = current position
parBestPos = parPos;
parBestFit = parFit;

fcount = 1; % Figure counter

if nFun > 3
    plotFun = 3;
else
    plotFun = nFun;
end

limits = zeros([plotFun 2]); % Initialise limits for graph output

% Print iteration, mean for every fitness function, numbber of PF particles
nonDomFitBar = mean(initPFDisp,1);
str = repmat('%3.4f ', 1, nFun);
fprintf(['Iteration %i: Mean PF f(x): ' str ' nPF: %i \n'], 0, nonDomFitBar, nNonDomParticles);

history(1,:) = [0, nonDomFitBar, nNonDomParticles];

%% Main PSO Loop
for it = 2:maxIt+1
    %% Select Global Best
    % Empty every timestep as it's size is altered below (repmat)
    GlobalBest = [];
    
    % Using two global bests, based on max and min particles dominated
    GlobalBest(1,:) = gbestminroulette(nonDomFit,parsDominated,nonDomPos);
    GlobalBest(2,:) = gbestmaxroulette(nonDomFit,parsDominated,nonDomPos);
    
    % Repeat Gbest to creat matrix with number of rows = nPop
    divisor = size(GlobalBest,1);
    repFactor = nPop/divisor;
    GlobalBest = repmat(GlobalBest,repFactor,1);
	
    %% Update currect position and enforce bounds/conditions
    % Update velocity > Update particle position
    parVel = w*parVel...
        + c1*rand(varSize).*(parBestPos - parPos)...
        + c2*rand(varSize).*(GlobalBest - parPos);
    
    parPos = parPos + parVel;
    
    % Mutation functions
    parPos(set2,:) = unimutation(parPos(set2,:),setSize,nVar,varMinMat,varMaxMat,mutProb);
    parPos(set3,:) = nonunimutation(parPos(set3,:),setSize,nVar,varMinMat(set3,:),varMaxMat(set3,:),mutProb,it-1,maxIt);
    
    % Ensure new particle position is within boundaries
    con1 = parPos > varMaxMat;
    con2 = parPos < varMinMat;
    
    % If not set particle velocity to zero and enforce bounds
    parVel(con1 | con2) = 0;
    parPos(con1) = varMaxMat(con1);
    parPos(con2) = varMinMat(con2);
    
    %% Only optimisation part dependent on problem
    
    % Calculate fitness functions
    [parFit,parPos] = costcaller(costFun,nPop,nFun,parPos,cond,varArray,options);
   
    %% Update particle best position
    % If new dominates old: new = best
    con1 = all(parFit <= parBestFit,2);
    t1 = nPopArray(con1);
    con1 = any(parFit(t1,:) < parBestFit(t1,:),2);
    t1 = t1(con1);
    parBestFit(t1,:) = parFit(t1,:);
	parBestPos(t1,:) = parPos(t1,:);
    
    % If old dominates new: no change
    con2 = all(parBestFit <= parFit,2);
    t2 = nPopArray(con2);
    con2 = any(parBestFit(t2,:) < parFit(t2,:),2);
    t2 = t2(con2);
    
    % If neither dominates, pick random
    con3 = ~ismember(nPopArray,[t1;t2]);
    t3 = nPopArray(con3);
    combineFit = [parFit(t3,:),parBestFit(t3,:)];
    combinePos = [parPos(t3,:),parBestPos(t3,:)];
    
    random = randi(2,length(t3),1);
    assignFit = random * nFun;
    assignPos = random * nVar;
    for i=1:length(assignFit)
        j = assignFit(i)-nFun+1:assignFit(i);
        parBestFit(t3(i),:) = combineFit(i,j);

        k = assignPos(i)-nVar+1:assignPos(i);
        parBestPos(t3(i),:) = combinePos(i,k);
    end
    
    %% Update Pareto Front
    % Combination of current particle positions and non-dominated positions
    combiPos = [parBestPos; nonDomPos];
    combiFit = [parBestFit; nonDomFit];
    combiUni = unique([combiPos,combiFit],'stable','rows');
    combiPos = combiUni(:,1:nVar);
    combiFit = combiUni(:,(1:nFun)+nVar);
    
    j=1;
    for i=1:size(combiFit,1)
        contender = combiFit(i,:);
        
        con1 = all(combiFit <= contender,2);
        con2 = any(combiFit(con1,:) < contender);
        
        if ~any(con2)
            con1 = all(contender <= combiFit,2);
            con2 = any(contender < combiFit(con1,:),2);
            parsDominated(j,:) = sum(con2);
            nonDomFit(j,:) = combiFit(i,:);
            nonDomPos(j,:) = combiPos(i,:);
            j=j+1;
        end
    end
    
    parsDominated(j:end,:) = [];
    nonDomFit(j:end,:) = [];
    nonDomPos(j:end,:) = [];
    
    %% Crowding Distance
    maxfi = max(nonDomFit,[],1);
    minfi = min(nonDomFit,[],1);
    
    nonDomFitNorm = (nonDomFit - minfi)./(maxfi - minfi);
    
    nonDomCrowd = crowder(nonDomFitNorm,maxfi,minfi,nFun);
    
    [nNonDomParticles,~] = size(nonDomCrowd);
    
    if nNonDomParticles > maxPF
        dim = (1:maxPF)';
        num = (1:size(nonDomCrowd,1))';
        csort = sortrows([num, nonDomCrowd],2,'descend');
        
        % Only need this if gbest selector needs it
        % nonDomCrowd = csort(dim,2);
    
        nonDomFit = nonDomFit(csort(dim,1),:);
        nonDomPos = nonDomPos(csort(dim,1),:);
        parsDominated = parsDominated(csort(dim,1));
        
        nNonDomParticles = maxPF;
        
    end
    
    %% Display
    parFitDisp = parFit;
    parFitDisp(:,inv) = inv(:,inv)./parFitDisp(:,inv);
    
    nonDomFitDisp = nonDomFit;
    nonDomFitDisp(:,inv) = inv(:,inv)./nonDomFitDisp(:,inv);
    
    nonDomFitBar = mean(nonDomFitDisp,1);
    fprintf(['Iteration %i: Mean PF f(x): ' str ' nPF: %i \n'], it-1, nonDomFitBar, nNonDomParticles);
    
    history(it,:) = [it-1, nonDomFitBar, nNonDomParticles];
    
    if mod(it-1,fi) == 0
        
        if options.Baseline
            
            baselineCost = options.Base.Results.Cost;
        else
            baselineCost = zeros(1,nFun);
        end
            
        % Have max limits slightly above that of PF so that points are not
        % lying on the edges
        maxf = max([nonDomFitDisp; baselineCost] ,[],1)*1.2;
        minf = min([nonDomFitDisp; baselineCost] ,[],1);
        figure(fcount)
        clf
        title(['Pareto Front (Iteration: ' num2str(it-1) ')'])
        set(gcf, 'Position', [0, 0, 1920, 1200])
        hold on
        xlabel('1/L')
        ylabel('Cd')
        
        % Max amount of cost functions that can be plotted is 3, so if
        % there are > 3, only plot first 3
        if nFun > 3
            maxf = maxf(1:3);
            minf = minf(1:3);
        end
            
        for i = 1:plotFun
            if inv(i) == 0
                limits(i,:) = [min(0,minf(i)), maxf(i)];
            else
                limits(i,:) = [minf(i), maxf(i)];
            end
        end
        
        xlim(limits(1,:))
        ylim(limits(2,:))
        
        if nFun == 2

            plot(baselineCost(:,1), baselineCost(:,2),'ko');
            % plot(parFitDisp(:,1), parFitDisp(:,2),'k*');
            plot(nonDomFitDisp(:,1), nonDomFitDisp(:,2),'kx');
        else
            plot3(baselineCost(:,1), baselineCost(:,2), baselineCost(:,3),'ko');
            % plot3(parFitDisp(:,1), parFitDisp(:,2), parFitDisp(:,3),'k*');
            plot3(nonDomFitDisp(:,1), nonDomFitDisp(:,2), nonDomFitDisp(:,3),'kx');
            zlim(limits(3,:))
            zlabel('M')
        end
        
        legend('Baseline', 'Optimal Design')
        hold off
        % Pause to display graph while simulation is running
        pause(0.00001)
        fcount = fcount + 1;
    end
    
    % Stops simulation if any Pareto Front values are not numeric 
    if any(any(isnan(nonDomFit)))
        
        error('Non-dominated particle(s) cost function value(s) non-numeric')
    end
end

%% Ordering for output
% Ordered with respect to fitness function 1, then 2, etc
num = (1:length(parsDominated))';

if nFun == 2
    
    fsort = sortrows([num, nonDomFit],2,'descend');
else
    fsort = sortrows([num, nonDomFit],[2 3],{'descend' 'descend'});
end

nonDomFit = nonDomFit(fsort(num,1),:);
nonDomPos = nonDomPos(fsort(num,1),:);

%% Different methods to find "ideal" Pareto Front point
% idealPoint = zeros(1,nFun);
% [fU,~] = max(nonDomFit);
% delta = (nonDomFit-idealPoint)./(fU-idealPoint);
% nonDomFitDist = sqrt(sum(delta.^2,2));
% [~,idxKP1] = min(nonDomFitDist);
% 
% particleviewer(nonDomPos(idxKP1,:));
% 
% [fU,~] = max(nonDomFit);
% [fL,~] = min(nonDomFit);
% delta = (nonDomFit-fL)./(fU-fL);
% nonDomFitDist = sqrt(sum(delta.^2,2));
% [~,idxKP2] = min(nonDomFitDist);
% 
% particleviewer(nonDomPos(idxKP2,:));