function area = panelarea(x,y,z)
%% Calculate panel area
% Trapezoidal area found by spliting into two triangles and summing, works
% for any type of triangle

ax = x(1:end-1,1:end-1);
bx = x(2:end,1:end-1);
cx = x(2:end,2:end);
dx = x(1:end-1,2:end);

ay = y(1:end-1,1:end-1);
by = y(2:end,1:end-1);
cy = y(2:end,2:end);
dy = y(1:end-1,2:end);

az = z(1:end-1,1:end-1);
bz = z(2:end,1:end-1);
cz = z(2:end,2:end);
dz = z(1:end-1,2:end);

% Areas based on xyz vectors so any triangle can be input
% Triangle one
base = sqrt((bx-ax).^2 + (by-ay).^2 + (bz-az).^2);
height = sqrt((dx-ax).^2 + (dy-ay).^2 + (dz-az).^2);

A1 = 0.5.*base.*height;

% Triangle two
base = sqrt((cx-dx).^2 + (cy-dy).^2 + (cz-dz).^2);
height = sqrt((cx-bx).^2 + (cy-by).^2 + (cz-bz).^2);

A2 = 0.5.*base.*height;

area = A1 + A2;