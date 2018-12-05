function run = flightstate(flow,i)
%% Flow conditions for specific flight state
run = flow;

Minf = flow.Minf(i);
alpha = flow.alpha(i) * pi/180;
a = flow.a(i);

Uinf = Minf * a;

planeAngles = [alpha 0 pi/2-alpha];
    
run.delq = flow.delq(i);
run.Machq = flow.Machq(i);

run.alpha = alpha;
run.Minf = Minf;

if ~isempty(flow.delta)
    run.delta = flow.delta(i) * pi/180;
end

run.U = Uinf*[cos(alpha) 0 sin(alpha)];
run.Uinf = Uinf;
run.Pinf = flow.Pinf(i);
run.rho = flow.rho(i);
run.Tinf = flow.Tinf(i);
run.mu = flow.mu(i);
run.a = flow.a(i);
run.Pr = flow.Pr(i);

% Angles between planes
run.planeAngles = planeAngles;