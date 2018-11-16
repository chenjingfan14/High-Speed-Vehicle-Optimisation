function flow = flowparameters(alpha,Mach,alt,delta)
%% Flow conditions 

if nargin > 0 && ~isempty(alpha)
else
    alpha = -4:4:24;
    alpha = 0;
end

if nargin > 1 && ~isempty(Mach)
else
    Mach = [4.63];
end

if nargin > 2 && ~isempty(alt)
else
    alt = [32850];
end

if nargin > 3 && ~isempty(delta)
else
    delta = [5];
end

alphaDim = numel(alpha);
MinfDim = numel(Mach);
altDim = numel(alt);
deltaDim = numel(delta);

parameterIndex = fullfact([alphaDim,MinfDim,altDim,deltaDim]);

[row,~] = size(parameterIndex);

gamma = 1.4;
R = 287;

% Freestream stagnation pressure ratio
Pinf_P0 = (2./((gamma+1)*(Mach.^2))).^(gamma/(gamma-1)) .* (((2*gamma*(Mach.^2))-(gamma-1))/(gamma+1)).^(1/(gamma-1));

for i = MinfDim:-1:1
    [matchdel(i),matchMach(i)] = matchingPoint(gamma,Pinf_P0(i));
end

for i = altDim:-1:1
    hvals = atmosphere(alt(i),0,0); % Find flow parameters
    
    Tinf(i) = hvals(1); % Freestream temperature
    rho(i) = hvals(2); % Freestream density
    Pinf(i) = hvals(7); % Freestream pressure
    a(i) = hvals(5); % Speed of sound
    mu(i) = hvals(8); % Dynamic viscosity
    kt(i) = hvals(10); % Thermal conductivity

    cp = R*gamma/(gamma-1);
    Pr(i) = mu*cp/kt;
end

flow.ParameterIndex = parameterIndex;
flow.alpha = alpha;
flow.Minf = Mach;
flow.delta = delta;
flow.gamma = gamma;
flow.Pinf = Pinf;
flow.Tinf = Tinf;
flow.R = R;
flow.Pr = Pr;
flow.rho = rho;
flow.mu = mu;
flow.a = a;
flow.delq = matchdel;
flow.Machq = matchMach;
flow.Dim = [row,1];
flow.Runs = row;