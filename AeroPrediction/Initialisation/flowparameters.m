function flow = flowparameters(AoA,M)
%% Flow conditions 

if nargin > 0 && ~isempty(AoA)
    alpha = AoA;
else
    alpha = -4:4:24;
end

if nargin > 1 && ~isempty(M)
    Minf = M;
else
    Minf = [4.63];
end

if ~isequal(size(alpha),size(Minf))
    MinfArray = reshape(Minf,1,[]);
    alphaArray = reshape(alpha,[],1);
    
    Minf = repmat(MinfArray,size(alphaArray));
    alpha = repmat(alphaArray,size(MinfArray));
end

[row,col] = size(alpha);

alt = 32850;
gamma = 1.4;
R = 287;

% Freestream stagnation pressure ratio
Pinf_P0 = (2./((gamma+1)*(Minf.^2))).^(gamma/(gamma-1)) .* (((2*gamma*(Minf.^2))-(gamma-1))/(gamma+1)).^(1/(gamma-1));

% Find Newtonian/PMe matching point
[dim1,dim2] = size(Minf);

for i = dim1:-1:1
    for j = dim2:-1:1
        [matchdel(i,j),matchMach(i,j)] = matchingPoint(gamma,Pinf_P0(i,j));
    end
end

hvals = atmosphere(alt,0,0); % Find flow parameters

a = hvals(5); % Speed of sound
mu = hvals(8); % Dynamic viscosity
kt = hvals(10); % Thermal conductivity
Uinf = Minf*a;

cp = R*gamma/(gamma-1);
Pr = mu*cp/kt;

flow.alpha = alpha;
flow.Minf = Minf;
flow.gamma = gamma;
flow.Uinf = Uinf;
flow.Pinf = hvals(7);
flow.Tinf = hvals(1);
flow.R = R;
flow.Pr = Pr;
flow.rho = hvals(2);
flow.mu = hvals(8);
flow.a = a;
flow.delq = matchdel;
flow.Machq = matchMach;
flow.runs = row*col;
flow.dim = [row,col];