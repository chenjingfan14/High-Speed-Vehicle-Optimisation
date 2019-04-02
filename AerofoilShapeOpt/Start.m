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
    
    % Begin Optimisation
    if method == "PSO"
        
        StartPSO()
        
    elseif method == "GA"
        
        StartGA()
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