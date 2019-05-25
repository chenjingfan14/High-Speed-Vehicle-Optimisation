function flow = flowparameters(inputs)
%% Flow conditions for all flight states
% Check if inputs have been provided
% If so, bypass. If not, give pre-determined values
if isfield(inputs,'Alpha')

    alpha = inputs.Alpha;
else
    alpha = (-4:4:24)';
end

if isfield(inputs,'Mach')

    Mach = inputs.Mach;
else
    Mach = [4.63]';
end

if isfield(inputs,'Altitude')

    alt = inputs.Altitude;
else
    alt = [32850]';
end

if isfield(inputs,'Control') 
    
    control = inputs.Control;
else
    control = false;
end

if isfield(inputs,'Delta')
    
    delta = inputs.Delta;
else
    delta = [0]';
end

%% Running Order
% Define order in which properties will be updated through the analysis

% Initial order ie. cycle all angle of attacks for same Mach, Altitude, 
% Control surface deflection, then increase Mach, repeat... etc.
flightStateOrder = ["Alpha","Mach","Altitude","Delta"];

% How each parameter will be plotted graphically
alphaPlot = ['Angle of Attack (' char(176) ')'];
deltaPlot = ['Control Surface Deflection (' char(176) ')'];
plotOrder = [alphaPlot, "Mach", "Altitude (ft)", deltaPlot];

% How many analysis points of each parameter
[alphaDim, alphaSets] = size(alpha);
MachDim = numel(Mach);
altDim = numel(alt);

if control % Same as above
    
    deltaDim = numel(delta);
else % Set to one so that it stays at the end and can be removed
    deltaDim = 1;
end

% Initial order (as above)
matIndex = [alphaDim,MachDim,altDim,deltaDim];

% Sort parameters from most analysis points to least. For example if only
% one Mach number is desired but 3 altitudes, better to output graphs in
% terms of "Mach number at these altitudes", rather than "Altitude 1 for
% Mach number", Altitude 2 for Mach number", etc.
[matIndex,order] = sort(matIndex,2,'descend');
flightStateOrder = flightStateOrder(order);
plotOrder = plotOrder(order);

%% Create flight state matrix
% Using new running order, create all flight state combinations and fill 
% with respective parameter value (Alpha, Mach, etc)

dim3Multi = 1;

if alphaSets > 1
    
    % Use for experimental comparisons where angle of attack may be
    % different across different Mach/altitude experiments
    flightStates(:,1) = alpha(:);
    flightStateIndex(:,1) = 1:(alphaDim * alphaSets);
    
    MachSets = repmat(reshape(Mach,1,[]),alphaDim,1);
    altSets = repmat(reshape(alt,1,[]),alphaDim,1);
    id = repmat(1:MachDim,alphaDim,1);
    
    flightStates(:,[2 3]) = [MachSets(:) altSets(:)];
    flightStateIndex(:,[2 3]) = [id(:) id(:)];
    
    [row,~] = size(flightStates);
    
    alphaCol = 1;
    MachCol = 2;
    altCol = 3;
    deltaCol = 4;
    
else
    % Number of variable flight conditions
    dim = numel(flightStateOrder);
    
    % Produces all index combinations of flight states
    flightStateIndex = fullfact(matIndex);
    
    % Initialise all flight state parameter combinations
    [row,col] = size(flightStateIndex);
    flightStates = zeros(row,col);
    
    for i = 1:dim
        switch flightStateOrder(i)
            case "Alpha"
                paramFill = alpha(flightStateIndex(:,i));
                alphaCol = i;
            case "Mach"
                paramFill = Mach(flightStateIndex(:,i));
                MachCol = i;
            case "Altitude"
                paramFill = alt(flightStateIndex(:,i));
                altCol = i;
            case "Delta"
                paramFill = delta(flightStateIndex(:,i));
                deltaCol = i;
        end
        
        flightStates(:,i) = paramFill;
        
        % Any dimensions >= 3 must be combined to produce 2D graphs (with
        % combined third dimension being number of 2D graphs for each set of
        % data)
        if i >= 3
            
            dim3Multi = dim3Multi * matIndex(i);
        end
    end   
end

% Remove control parameters if no control specified
if ~control
    flightStateOrder(deltaCol) = [];
    
    if alphaSets == 1
        
        flightStateIndex(:,deltaCol) = [];
        flightStates(:,deltaCol) = [];
    end
    
    deltaCol = [];
end

% New output dimensions
dim = [matIndex(1:2),dim3Multi];

%% Mach/Altitude dependent properties
gamma = 1.4;
R = 287;

% Freestream stagnation pressure ratio
Pinf_P0 = (2./((gamma+1)*(Mach.^2))).^(gamma/(gamma-1)) .* (((2*gamma*(Mach.^2))-(gamma-1))/(gamma+1)).^(1/(gamma-1));

% Matching points for Newtonian + Prandtl-Meyer method
% CHECK: reasonable results for all Mach numbers?
for i = MachDim:-1:1
    [matchdel(i,:),matchMach(i,:)] = matchingPoint(gamma,Pinf_P0(i));
end

for i = altDim:-1:1
    hvals = atmosphere(alt(i),0,0); % Find flow parameters
    
    Tinf(i,:) = hvals(1); % Freestream temperature
    rho(i,:) = hvals(2); % Freestream density
    Pinf(i,:) = hvals(7); % Freestream pressure
    a(i,:) = hvals(5); % Speed of sound
    mu(i,:) = hvals(8); % Dynamic viscosity
    kt(i,:) = hvals(10); % Thermal conductivity

end

cp = R*gamma/(gamma-1);
Pr = mu.*cp./kt; % Prandtl number

%% Save flow parameters for all flight conditions
flow.FlightStateOrder = flightStateOrder;
flow.FlightStates = flightStates;
flow.PlotOrder = plotOrder;
flow.Runs = row;
flow.Dim = dim;

% Save flight conditions with respect to what they depend on
flow.alpha = alpha(flightStateIndex(:,alphaCol));
flow.Minf = Mach(flightStateIndex(:,MachCol));
flow.Altitude = alt(flightStateIndex(:,altCol));
flow.delta = delta(flightStateIndex(:,deltaCol));
flow.gamma = gamma;
flow.Pinf = Pinf(flightStateIndex(:,altCol));
flow.Tinf = Tinf(flightStateIndex(:,altCol));
flow.R = R;
flow.Pr = Pr(flightStateIndex(:,altCol));
flow.rho = rho(flightStateIndex(:,altCol));
flow.mu = mu(flightStateIndex(:,altCol));
flow.a = a(flightStateIndex(:,altCol));
flow.Uinf = flow.Minf .* flow.a;
flow.q = 0.5 * flow.rho .* flow.Uinf.^2;
flow.delq = matchdel(flightStateIndex(:,MachCol));
flow.Machq = matchMach(flightStateIndex(:,MachCol));