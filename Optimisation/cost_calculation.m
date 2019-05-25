function cost = cost_calculation(results,parameters,options)

nFun = options.CostFunctions;
wing = options.Wing;
baseline = options.Baseline;

if wing
    
    % lift = Lbar/wingspan;
    % moment = Mbar/wingspan;
    
    % Constraint section. Apply penalties if desired values are too
    % high/low
    
    if baseline
        
        base = options.Base;
        baseResults = base.Results;
        baseParameters = base.Parameters;
        
        constrain = [parameters.Aref/2, parameters.Wingspan, ...
            min(parameters.Sweep), max(parameters.Sweep), results.Mbar, ...
            results.copBar, results.Lbar, max(parameters.ThicknessND), max(parameters.Thickness)];
        
        minVal = [baseParameters.Aref*0.5, 0, 0, 0, 0, 0, baseResults.Lbar, 0, 0];
        
        maxVal = [baseParameters.Aref*2, baseParameters.Wingspan, 80, 80, baseResults.Mbar, ...
            baseResults.copBar, baseResults.Lbar*10, max(baseParameters.ThicknessND), max(baseParameters.Thickness)];
    else
        constrain = [results.Mbar,results.copBar];
        minVal = [0,0];
        maxVal = [inf,inf];
    end
    
    penalty = violation(constrain,minVal,maxVal);
    
    % Check magnitude of penalty
%     if ~any(penalty)
%         
%         penalty
%     end
    
    %% Costs
    % Varies automatically depending on number of cost functions defined in
    % simOptions
    
    if nFun == 2
        
        % Multi-objective
%         cost = [1/results.Lbar, results.Cdbar] + penalty;
        cost = [-results.Clbar, results.Cdbar] + penalty;
        
    elseif nFun == 1
        
        % Single objective
        cost = results.Cdbar/results.Clbar + penalty;
        
    else
        error(['No option for ' num2str(nFun) ' cost functions'])
    end
else
    cost = [];
end