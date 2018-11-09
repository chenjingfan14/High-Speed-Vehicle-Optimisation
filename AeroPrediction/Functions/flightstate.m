function run = flightstate(flow,i)
%% Flow conditions 

run = flow;

alpha = flow.alpha(i) * pi/180;
Uinf = flow.Uinf(i);

planeAngles = [alpha 0 pi/2-alpha];
    
run.delq = flow.delq(i);
run.Machq = flow.Machq(i);

run.alpha = alpha;
run.Minf = flow.Minf(i);
run.U = Uinf*[cos(alpha) 0 sin(alpha)];
run.Uinf = Uinf;

% Angles between planes
run.planeAngles = planeAngles;