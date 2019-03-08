function partStruct = normals(points)
%% Additional panel specific properties
% Calculate centre, outward facing normals, area

% Initialise and calculate centre points
a = points(1:end-1,1:end-1,:);
b = points(2:end,1:end-1,:);
c = points(2:end,2:end,:);
d = points(1:end-1,2:end,:);

centre = (a + b + c + d)/4;

radialLocation = atan2d(centre(:,:,2) - points(1,1,2), centre(:,:,3) - points(1,1,3));

%% Calculate normals of each panel
ac = c - a;
db = b - d;

% Normal vector and its magnitude
norm = crossmat(db, ac);
magNorm = magmat(norm);

unitNorm = norm./magNorm;

% Some normal magnitudes are zero which leads to NaN. Replace such
% values with zero
con = isnan(unitNorm);
unitNorm(con) = 0;

nx = unitNorm(:,:,1);
nyz = (unitNorm(:,:,2).^2 + unitNorm(:,:,3).^2).^0.5;
halfAngle = atan2(-nx, nyz) * 180/pi;

flow = nx < 0;

del = round(halfAngle,10);

con1 = del > 90;
con2 = del < -90;
del(con1) = 180 - del(con1);
del(con2) = -180 - del(con2);
del = del*pi/180;

%% Calculate average deltas across each panel

% Might not be needed if not using Pollock's method for viscous
% particle tracking
xDiff = diff(points,1,1);
yDiff = diff(points,1,2);

deltax = (xDiff(:,:,1).^2 + xDiff(:,:,2).^2 + xDiff(:,:,3).^2).^0.5;
deltay = (yDiff(:,:,1).^2 + yDiff(:,:,2).^2 + yDiff(:,:,3).^2).^0.5;

%% Calculate panel areas
% Areas based on xyz vectors so any triangle can be input
% Triangle one
base = sqrt((b(:,:,1) - a(:,:,1)).^2 + (b(:,:,2) - a(:,:,2)).^2 + (b(:,:,3) - a(:,:,3)).^2);
height = sqrt((d(:,:,1) - a(:,:,1)).^2 + (d(:,:,2) - a(:,:,2)).^2 + (d(:,:,3) - a(:,:,3)).^2);

A1 = 0.5.*base.*height;

% Triangle two
base = sqrt((c(:,:,1) - d(:,:,1)).^2 + (c(:,:,2) - d(:,:,2)).^2 + (c(:,:,3) - d(:,:,3)).^2);
height = sqrt((c(:,:,1) - b(:,:,1)).^2 + (c(:,:,2) - b(:,:,2)).^2 + (c(:,:,3) - b(:,:,3)).^2);

A2 = 0.5.*base.*height;

area = A1 + A2;

partStruct.Points = points;
partStruct.a = a;
partStruct.b = b;
partStruct.c = c;
partStruct.d = d;
partStruct.centre = centre;
partStruct.radialLocation = radialLocation;
partStruct.norm = norm;
partStruct.unitNorm = unitNorm;
partStruct.deltax = deltax;
partStruct.deltay = deltay;
partStruct.area = area;
partStruct.del = del;
partStruct.flow = flow;

end