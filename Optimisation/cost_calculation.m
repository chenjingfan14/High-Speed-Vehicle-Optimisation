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
        baseResults = options.Base.Results;
        
        constrain = [min(parameters.Sweep),results.Mbar,results.copBar,results.Lbar];
        minVal = [0,0,0,baseResults.Lbar];
        maxVal = [80,baseResults.Mbar,baseResults.copBar,inf];
    else
        constrain = [results.Mbar,results.copBar];
        minVal = [0,0];
        maxVal = [inf,0.6];
    end
    
    penalty = violation(constrain,minVal,maxVal);
    
    % Check magnitude of penalty
    %     if any(penalty)
    %         penalty
    %     end
    
    %% Costs
    % Varies automatically depending on number of cost functions defined in
    % simOptions
    
    if nFun == 2
        
        % Multi-objective
        cost = [1/results.Lbar,results.Cdbar] + penalty;
        
    elseif nFun == 1
        
        % Single objective
        cost = results.Cdbar/results.Clbar + penalty;
        
    else
        error(['No option for ' num2str(nFun) ' cost functions'])
    end
    
    % If any cost less than zero, particle swarm will see it as optimal,
    % whereas none of these aerodynamic values should be less than zero
    infCon = cost < 0;
    
    cost(infCon) = inf;
else
    cost = [];
end