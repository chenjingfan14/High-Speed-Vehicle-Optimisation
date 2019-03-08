function dt = estimateexit(xp1,yp1,Vp,xc,yc,xPrevCross,yPrevCross,mini,maxi)

x1 = xc(1,:);
x2 = xc(2,:);
x3 = xc(3,:);
x4 = xc(4,:);

y1 = yc(1,:);
y2 = yc(2,:);
y3 = yc(3,:);
y4 = yc(4,:);

dtLarge = 10;

xp2 = xp1 + Vp(1,:)*dtLarge;
yp2 = yp1 + Vp(2,:)*dtLarge;

dimx = size(xPrevCross);
dimy = size(yPrevCross);

xCross = nan(dimx);
yCross = nan(dimy);

%% Has particle crossed between panel points 1 & 2
% Intersection point of two infinite lines
Px1 = ((xp1.*yp2 - yp1.*xp2).*(x1-x2) - (xp1-xp2).*(x1.*y2 - y1.*x2))./...
    ((xp1-xp2).*(y1-y2) - (yp1-yp2).*(x1-x2));

Py1 = ((xp1.*yp2 - yp1.*xp2).*(y1-y2) - (yp1-yp2).*(x1.*y2 - y1.*x2))./...
    ((xp1-xp2).*(y1-y2) - (yp1-yp2).*(x1-x2));

% Defines whether intersection point is between two given line segments
t1 = ((xp1-x1).*(y1-y2) - (yp1-y1).*(x1-x2))./...
    ((xp1-xp2).*(y1-y2) - (yp1-yp2).*(x1-x2));

u1 = -((xp1-xp2).*(yp1-y1) - (yp1-yp2).*(xp1-x1))./...
    ((xp1-xp2).*(y1-y2) - (yp1-yp2).*(x1-x2));

% Is intersection is between both line segments
ydec = t1 >= mini & t1 <= maxi & u1 >= mini & u1 <= maxi;

% Which panel line is being crossed
yCross(:,ydec) = [x1(ydec);y1(ydec);x2(ydec);y2(ydec)];

% Particle cannot cross back over line it has just crossed, must be
% different otherwise no intersection
ydec = ydec & any(yCross ~= yPrevCross,1);

%% Has particle crossed between panel points 2 & 3
Px2 = ((xp1.*yp2 - yp1.*xp2).*(x2-x3) - (xp1-xp2).*(x2.*y3 - y2.*x3))./...
    ((xp1-xp2).*(y2-y3) - (yp1-yp2).*(x2-x3));

Py2 = ((xp1.*yp2 - yp1.*xp2).*(y2-y3) - (yp1-yp2).*(x2.*y3 - y2.*x3))./...
    ((xp1-xp2).*(y2-y3) - (yp1-yp2).*(x2-x3));

t2 = ((xp1-x2).*(y2-y3) - (yp1-y2).*(x2-x3))./...
    ((xp1-xp2).*(y2-y3) - (yp1-yp2).*(x2-x3));

u2 = -((xp1-xp2).*(yp1-y2) - (yp1-yp2).*(xp1-x2))./...
    ((xp1-xp2).*(y2-y3) - (yp1-yp2).*(x2-x3));

xinc = t2 >= mini & t2 <= maxi & u2 >= mini & u2 <= maxi;
xCross(:,xinc) = [x2(xinc);y2(xinc);x3(xinc);y3(xinc)];
xinc = xinc & any(xCross ~= xPrevCross,1);

%% Has particle crossed between panel points 3 & 4
Px3 = ((xp1.*yp2 - yp1.*xp2).*(x3-x4) - (xp1-xp2).*(x3.*y4 - y3.*x4))./...
    ((xp1-xp2).*(y3-y4) - (yp1-yp2).*(x3-x4));

Py3 = ((xp1.*yp2 - yp1.*xp2).*(y3-y4) - (yp1-yp2).*(x3.*y4 - y3.*x4))./...
    ((xp1-xp2).*(y3-y4) - (yp1-yp2).*(x3-x4));

t3 = ((xp1-x3).*(y3-y4) - (yp1-y3).*(x3-x4))./...
    ((xp1-xp2).*(y3-y4) - (yp1-yp2).*(x3-x4));

u3 = -((xp1-xp2).*(yp1-y3) - (yp1-yp2).*(xp1-x3))./...
    ((xp1-xp2).*(y3-y4) - (yp1-yp2).*(x3-x4));

yinc = t3 >= mini & t3 <= maxi & u3 >= mini & u3 <= maxi;
yCross(:,yinc) = [x4(yinc);y4(yinc);x3(yinc);y3(yinc)];
yinc = yinc & any(yCross ~= yPrevCross,1);

xpt(ydec) = Px1(ydec);
xpt(xinc) = Px2(xinc);
xpt(yinc) = Px3(yinc);

ypt(ydec) = Py1(ydec);
ypt(xinc) = Py2(xinc);
ypt(yinc) = Py3(yinc);

p1 = (xp1.^2 + yp1.^2).^0.5;
pt = (xpt.^2 + ypt.^2).^0.5;

Vmag = (Vp(1,:).^2 + Vp(2,:).^2).^0.5;

dt = (abs(pt - p1))./Vmag;

dt = dt/2;

