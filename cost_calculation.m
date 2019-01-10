function cost = cost_calculation(results,options)

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
        
        constrain = [results.Mbar,results.copBar,results.Lbar];
        minVal = [0,0,baseResults.Lbar];
        maxVal = [baseResults.Mbar,baseResults.copBar,inf];
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
    
    cost = [1/results.Lbar,results.Cdbar] + penalty;
    
    % If any cost less than zero, particle swarm will see it as optimal,
    % whereas none of these aerodynamic values should be less than zero
    infCon = cost < 0;
    
    cost(infCon) = inf;
    
else
    cost = [];
end