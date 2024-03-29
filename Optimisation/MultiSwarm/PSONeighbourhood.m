function [GlobalBestFit,GlobalBestPos,history] = PSONeighbourhood(cond,costFun,varArray,varMin,varMax,nVar,nPop,maxIt,maxStall,mutProb,w,wmax,wmin,c1,c2,nFun,fi,options)

Ni = floor(0.25*nPop);
N = Ni;

tolerance = options.Tolerance;

% Creating 3 sets within population, defines mutation method on set
% However all particles still communicate with each other regardless of
% sets
setSize = nPop/3;
set1 = 1:setSize; % No mutation
set2 = set1 + setSize; % Uniform mutation
set3 = set2 + setSize; % Non-uniform mutation

%% Swarm

% Maximum velocity particle can have per iteration
varMinMat = repmat(varMin',nPop,1);
varMaxMat = repmat(varMax',nPop,1);
Vmax = 0.5*(varMaxMat-varMinMat);

varSize = [nPop nVar];

history = zeros(maxIt, 3);

%% Initialise

% Array to hold best cost value on each iteration
GlobalBestFit = zeros(maxIt+1,1);
GlobalBestPos = zeros(maxIt+1,nVar);

% Create population array
parPos = unifrnd(varMinMat,varMaxMat,varSize);
parVel = zeros(varSize);

% Calculate initial fitness functions
[parFit,parPos] = costcaller(costFun,nPop,nFun,parPos,cond,varArray,options);

parBestPos = parPos;
parBestFit = parFit;
[GlobalBestFit(1),bestID] = min(parBestFit);
GlobalBestPos(1,:) = parBestPos(bestID,:);

fcount = 1; % Figure counter

inv = options.Inv;
neg = options.Neg;

parFitDisp = parFit;
parFitDisp(:,inv) = inv(:,inv)./parFitDisp(:,inv);
parFitDisp(:,neg) = -parFitDisp(:,neg);

GlobalBestFitDisp = GlobalBestFit;
GlobalBestFitDisp(:,inv) = inv(:,inv)./GlobalBestFitDisp(:,inv);
GlobalBestFitDisp(:,neg) = -GlobalBestFitDisp(:,neg);

% Print iteration, mean for every fitness function, success rate
parFitBar = mean(parFitDisp(parFitDisp < inf),1);
fprintf('Iteration 0: Global Best: %3.4f Mean: %3.4f Stall: 0 \n', GlobalBestFitDisp(1), parFitBar);

history(1,:) = [0, GlobalBestFitDisp(1), 0];

%% Main PSO Loop
stall = 0;
wcounter = 0;

for it = 2:maxIt+1
    
    % Create random neighbourhoods based on variable neighbourhood size
    hood = zeros(nPop, N);
    hood(:,1) = 1:nPop;
    
    for i=1:nPop
        
        neighbours = randperm(nPop-1,N-1);
        iShift = neighbours >= i;
        neighbours(iShift) = neighbours(iShift) + 1;
        hood(i,2:end) = neighbours;
        
    end
    
    % Find best fitness value in neighbourhood
    hoodFit = parBestFit(hood);
    [hoodBestFit,ind] = min(hoodFit,[],2);
    vec = 1:size(hoodBestFit,1);
    idx = vec' + (ind-1)*length(vec);
    ind = hood(idx);
    hoodBestPos = parBestPos(ind,:);
    
    %% Update particles and enforce conditions
    parVel = w*parVel...
        + c1*rand(varSize).*(parBestPos - parPos)...
        + c2*rand(varSize).*(hoodBestPos - parPos);
 
    con1 = parVel > Vmax;
    con2 = parVel < -Vmax;
    
    parVel(con1) = Vmax(con1);
    parVel(con2) = -Vmax(con2);
    
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
    
    % Calculate fitness functions
    [parFit,parPos] = costcaller(costFun,nPop,nFun,parPos,cond,varArray,options);
    
    %% Update best particle/global fitness values and corresponding positions
    ind = parFit < parBestFit;
    parBestFit(ind) = parFit(ind);
    parBestPos(ind,:) = parPos(ind,:);
    
    [BestFit,id] = min(parBestFit);
    
    % If new best value, update global best, set stall counter to zero, set
    % neighbourhood size to minimum, reduce inertia counter
    % Otherwise, global best = previous value, increase stall counter, if
    % max stall value reached then end particle swarm, increase
    % neighbourhood size incrementally with initial value, increase inertia
    % counter
    if BestFit < GlobalBestFit(it-1)
        GlobalBestFit(it) = BestFit;
        GlobalBestPos(it,:) = parBestPos(id,:);
        stall = 0; 
        N = Ni;
        wcounter = max(0, wcounter-1);
    else
        GlobalBestFit(it) = GlobalBestFit(it-1);
        GlobalBestPos(it,:) = GlobalBestPos(it-1,:);
        stall = stall+1;
        
        N = min(N+Ni, nPop);
        wcounter = wcounter+1;
    end
    
    % Alter inertia based on above counter
    if wcounter < 2
        w=min(w*2,wmax);
    elseif wcounter > 5
        w=max(w*0.5,wmin);
    end
    
    %% Display
    parFitDisp = parFit;
    parFitDisp(:,inv) = inv(:,inv)./parFitDisp(:,inv);
    parFitDisp(:,neg) = -parFitDisp(:,neg);
    
    GlobalBestFitDisp = GlobalBestFit;
    GlobalBestFitDisp(:,inv) = inv(:,inv)./GlobalBestFitDisp(:,inv);
    GlobalBestFitDisp(:,neg) = -GlobalBestFitDisp(:,neg);
    
    parFitBar = mean(parFitDisp(parFitDisp < inf),1);
    fprintf('Iteration %i: Global Best: %3.4f Mean: %3.4f Stall: %i \n', it-1, GlobalBestFitDisp(it), parFitBar, stall);
    
    history(it,:) = [it-1, GlobalBestFitDisp(it), stall];
    
    convergence = stall == maxStall || GlobalBestFit(it) < tolerance;
    
    if mod(it-1,fi) == 0 || convergence
        
        singlefunplot(it-1,GlobalBestFitDisp,options,fcount)
        
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