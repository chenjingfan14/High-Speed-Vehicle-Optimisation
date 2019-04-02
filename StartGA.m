%% Start Genetic Algorithm

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

[GlobalBestFit,GlobalBestPos,history] = GeneticAlgorithm(cond,costFun,varArray,varMin,varMax,nVar,nPop,maxIt,maxStall,Pc,Pm,Er,nFun,inv,fi,options);

time = toc;