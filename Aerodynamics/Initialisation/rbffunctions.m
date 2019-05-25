function RBFfun = rbffunctions(method)

switch method
    
    case {"TPS","Thin Plate Spline"}
        
        RBFfun = @(X)X.^2 .* log(X);
        
    case "Wendlands"
        
        RBFfun = @(X)(1 - X).^2;
        
    case "Euclids Hat"
        
        % r = 200
        RBFfun = @(X)pi * (((1/12) .* X.^3) - 200.^2 .* X + (4/3 .* 200.^3));
        
    otherwise
        
        warning('RBF function set as default: Volume Spline')
        RBFfun = @(X)X;
end
        
        