function run = flightstate(flow,i)
%% Flow conditions 

alpha = flow.alpha(i) * pi/180;
Uinf = flow.Uinf(i);

planeAngles = [alpha 0 pi/2-alpha];
    
run.delq = flow.delq(i);
run.Machq = flow.Machq(i);

run.alpha = alpha;
run.Minf = flow.Minf(i);
run.gamma = flow.gamma;
run.U = Uinf*[cos(alpha) 0 sin(alpha)];
run.Uinf = Uinf;
run.Pinf = flow.Pinf;
run.rho = flow.rho;
run.a = flow.a;

% Angles between planes
run.planeAngles = planeAngles;