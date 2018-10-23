clear all
clc

l = 1.764*0.174;
r1 = 0.175*0.174;

xPanels = 200;
x = 0:l/(xPanels):l;

yzPanels = 90;
theta = 0:pi/yzPanels:pi;

breakpoint = 0.15*0.174;

con = x/l <= breakpoint/l;
con2 = x/l >= breakpoint/l;
xp = x(con);
xp2 = x(con2);

if xp(end) ~= breakpoint
    xp = [xp,breakpoint];
end
if xp2(1) ~= breakpoint
    xp2 = [breakpoint,xp2];
end

part = ((r1)^2 - (xp-r1).^2).^0.5;
part2 = (part(end) + 0.203*((xp2)-xp(end)));

y = round(part'.*sin(theta),4);
z = round(part'.*cos(theta),4);

y2 = round(part2'.*sin(theta),4);
z2 = round(part2'.*cos(theta),4);

y2(1,:) = y(end,:);
z2(1,:) = z(end,:);

config(1).Name = "test";
config(1).x = xp';
config(1).y = y;
config(1).z = z;

config(2).Name = "test";
config(2).x = xp2';
config(2).y = y2;
config(2).z = z2;

config = xyztopoints(config);

Properties(1).Conical = 1;
Properties(2).Conical = 1;
Properties(1).Points = config(1);
Properties(2).Points = config(2);
cellprop{1} = Properties(1);
cellprop{2} = Properties(2);

flow = flowparameters;
load('thetaBetaCurves.mat');
load('PrandtlMeyerExpansion.mat');

Aref = pi*max(part2)^2;

costFun = @aeroprediction;

plotter(config)

result = feval(costFun,cellprop,flow,Aref,0,thetaBetaM,maxThetaBetaM,PrandtlMeyer);
