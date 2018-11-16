function run = flightstate(flow,i)
%% Flow conditions 

run = flow;

index = flow.ParameterIndex;

alphaDep = index(i,1);
MachDep = index(i,2);
altDep = index(i,3);
deltaDep = index(i,4);

Minf = flow.Minf(MachDep);
alpha = flow.alpha(alphaDep) * pi/180;
a = flow.a(altDep);

Uinf = Minf * a;

planeAngles = [alpha 0 pi/2-alpha];
    
run.delq = flow.delq(MachDep);
run.Machq = flow.Machq(MachDep);

run.alpha = alpha;
run.Minf = Minf;
run.delta = flow.delta(deltaDep) * pi/180;
run.U = Uinf*[cos(alpha) 0 sin(alpha)];
run.Uinf = Uinf;
run.Pinf = flow.Pinf(altDep);
run.rho = flow.rho(altDep);
run.Tinf = flow.Tinf(altDep);
run.mu = flow.mu(altDep);
run.a = flow.a(altDep);
run.Pr = flow.Pr(altDep);

% Angles between planes
run.planeAngles = planeAngles;