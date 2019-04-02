%% Start Optimisation

% Initialise problem specific parameters, adds them to overall workspace
problemdefinition();

method = "PSO";
foils = getaerofoilsecdata;

tic
    
% Begin Optimisation
if method == "PSO"
    
    StartPSO()
    
elseif method == "GA"
    
    StartGA()
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