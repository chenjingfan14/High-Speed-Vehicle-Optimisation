clear
close all

% Use for incoming flow
Vinf = [10,0,0];

% Use for interpolating 4 corner plate velocities
% NOTE: Corners here defined as anti-clockwise from bottom left corner
Vc1 = [100,20];
Vc2 = [100,30];
Vc3 = [100,40];
Vc4 = [100,50];

%% Creating uniform/randomly perturbed panels
% Replace with size(points.x)
x = ones(10,50); % x,y points (x-1,y-1 panels)

% Creating x,y point matrices from 0 to [x y]-1
% Number of points
[dim1,dim2] = size(x);

% Number of panels
dim3 = dim1-1;
dim4 = dim2-1;

x = repmat((0:dim1-1)',1,dim2);
y = repmat(0:dim2-1,dim1,1);
z = zeros(dim1,dim2);

% Random perturbation to make uniform grid chaotic, only inner points
randx = rand(dim3-1,dim4-1)*0.5 - 0.25;
randy = rand(dim3-1,dim4-1)*0.5 - 0.25;

% randx = rand(dim1-2,dim2-2) - 0.5;
% randy = rand(dim1-2,dim2-2) - 0.5;

x(2:end-1,2:end-1) = x(2:end-1,2:end-1) + randx;
y(2:end-1,2:end-1) = y(2:end-1,2:end-1) + randy;

% Create points matrix to feed into normals function
points(:,:,1) = x;
points(:,:,2) = y;
points(:,:,3) = z;

points = normals(points);
plotter(points)

norm = points.unitNorm;

nx = norm(:,:,1);
ny = norm(:,:,2);
nz = norm(:,:,3);

a = points.a;
b = points.b;
c = points.c;
d = points.d;

ax = a(:,:,1);
bx = b(:,:,1);
cx = c(:,:,1);
dx = d(:,:,1);

ay = a(:,:,2);
by = b(:,:,2);
cy = c(:,:,2);
dy = d(:,:,2);

% Edge mid points
chordx = (x(:,2:end) + x(:,1:end-1))/2;
chordy = (y(:,2:end) + y(:,1:end-1))/2;

spanx = (x(2:end,:) + x(1:end-1,:))/2;
spany = (y(2:end,:) + y(1:end-1,:))/2;

% Initial injection points
xinject = [chordx(1,:)]; %spanx(:,1)'];
yinject = [chordy(1,:)]; %spany(:,1)'];

%%
rows = ones(1,dim4);
cols = 1:dim4;

rows = [rows]; %1:dim3];
cols = [cols]; %ones(1,dim3)];

nRows = dim3*100;

% Initial matrix ID, first leading edge panels
ID = rows + (cols-1)*dim3;
[~,tdim] = size(xinject);
tcols = 1:tdim;
timeID = 2 + (tcols-1)*nRows;

dummy = nan(nRows,tdim);

xp = dummy;
yp = dummy;

x1 = ax(ID);
x2 = bx(ID);
x3 = cx(ID);
x4 = dx(ID);

y1 = ay(ID);
y2 = by(ID);
y3 = cy(ID);
y4 = dy(ID);

c1 = [x1; y1];
c2 = [x2; y2];
c3 = [x3; y3];
c4 = [x4; y4];

%% Velocity definition

% Use for interpolating velocity from centre of panel to edge midpoints for
% uniform grids (potentially unit square?)
% [xV,yV,Vx1,Vx2,Vx3,Vx4,Vy1,Vy2,Vy3,Vy4] = centretomidvelocity(points,Vinf);

% Use for interpolating surface corner to velocities to panel corners for
% uniform grids. Both methods give basically same result

%% Bilinear interpolation, full form, no need to translate into unit square
mat1 = [1, x(1), y(1), x(1)*y(1);
    1, x(1), y(end), x(1)*y(end);
    1, x(end), y(1), x(end)*y(1);
    1, x(end), y(end), x(end)*y(end)];

xtemp = reshape(x,1,[]);
ytemp = reshape(y,1,[]);

mat2 = [ones(1,dim1*dim2); xtemp; ytemp; xtemp.*ytemp];

b = inv(mat1)'*mat2;

Vx = b(1,:)*Vc1(1) + b(2,:).*Vc4(1) + b(3,:).*Vc2(1) + b(4,:).*Vc3(1);
Vy = b(1,:)*Vc1(2) + b(2,:).*Vc4(2) + b(3,:).*Vc2(2) + b(4,:).*Vc3(2);

% Vx = b(1,:)*Vc1(1) + b(2,:).*Vc2(1) + b(3,:).*Vc3(1) + b(4,:).*Vc4(1);
% Vy = b(1,:)*Vc1(2) + b(2,:).*Vc2(2) + b(3,:).*Vc3(2) + b(4,:).*Vc4(2);

Vx = reshape(Vx,dim1,[]);
Vy = reshape(Vy,dim1,[]);

%% Translate to panel corner velocities
Vx1 = Vx(1:end-1,1:end-1);
Vx2 = Vx(2:end,1:end-1);
Vx3 = Vx(2:end,2:end);
Vx4 = Vx(1:end-1,2:end);

Vy1 = Vy(1:end-1,1:end-1);
Vy2 = Vy(2:end,1:end-1);
Vy3 = Vy(2:end,2:end);
Vy4 = Vy(1:end-1,2:end);

%% Initialise integration

stopCon = false(1,tdim);

timestep = 0.005;

t = 1;

V1 = [Vx1(ID); Vy1(ID)];
V2 = [Vx2(ID); Vy2(ID)];
V3 = [Vx3(ID); Vy3(ID)];
V4 = [Vx4(ID); Vy4(ID)];

xp(1,:) = xinject;
yp(1,:) = yinject;

Vp1 = V1+(yinject-c1).*((V4-V1)./(y4-y1));

mini = -0.00001;
maxi = 1.00001;

xp1 = xp(1,:);
yp1 = yp(1,:);
xp2 = xp(1,:) + 1e-16;
yp2 = yp(1,:) + 1e-16;

xCross = nan(4,tdim);
yCross = nan(4,tdim);
xPrevCross = nan(4,tdim);
yPrevCross = nan(4,tdim);

%% 2D
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

%% Updates for next timestep
% Line crossed in current timestep becomes previous crossing point for
% future timesteps
xPrevCross(:,xinc) = xCross(:,xinc);
yPrevCross(:,ydec) = yCross(:,ydec);
yPrevCross(:,yinc) = yCross(:,yinc);

realtime = zeros(1,tdim);
interpmodel = 'barycentricmodel';
particleArray = 1:tdim;
integrate = true;

if integrate
    while any(~stopCon)
        
        pt = [xp1;yp1];
        
        % RK4 integration
        k1 = timestep*feval(interpmodel, [pt; Vp1] ,timestep,x1,x2,x3,x4,y1,y2,y3,y4,V1,V2,V3,V4);
        k2 = timestep*feval(interpmodel, [pt; Vp1] + k1/2 ,timestep,x1,x2,x3,x4,y1,y2,y3,y4,V1,V2,V3,V4);
        k3 = timestep*feval(interpmodel, [pt; Vp1] + k2/2 ,timestep,x1,x2,x3,x4,y1,y2,y3,y4,V1,V2,V3,V4);
        k4 = timestep*feval(interpmodel, [pt; Vp1] + k3 ,timestep,x1,x2,x3,x4,y1,y2,y3,y4,V1,V2,V3,V4);
        
        delta = (k1 + 2*k2 + 2*k3 + k4)/6;
        
        deltax = delta(1,:);
        deltay = delta(2,:);
        
        Vp2 = Vp1 + delta(3:4,:);
        
        xpt = xp1 + deltax;
        ypt = yp1 + deltay;
        
        xp2 = xpt;
        yp2 = ypt;
        
        %% 2D
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
        
        %% Updates for next timestep
        % Line crossed in current timestep becomes previous crossing point for
        % future timesteps
        xPrevCross(:,xinc) = xCross(:,xinc);
        yPrevCross(:,ydec) = yCross(:,ydec);
        yPrevCross(:,yinc) = yCross(:,yinc);
        
        % Updates panel that particle is within
        rows(xinc) = rows(xinc) + 1;
        cols(yinc) = cols(yinc) + 1;
        cols(ydec) = cols(ydec) - 1;
        
        % Starting point for next timestep is intersection point. Removes error
        % of particle crossing current panel AND a t+1 panel line which will
        % not be picked up by simulation. However time and velocity must now be
        % updated for such situations
        xpt(ydec) = Px1(ydec);
        xpt(xinc) = Px2(xinc);
        xpt(yinc) = Px3(yinc);
        
        ypt(ydec) = Py1(ydec);
        ypt(xinc) = Py2(xinc);
        ypt(yinc) = Py3(yinc);
        
        % Streamline points updated
        xp(timeID) = xpt;
        yp(timeID) = ypt;
        
        % Plotting every timestep
%         figure(1)
%         hold on
%         for i=1:sum(~stopCon)
%             dim = timeID(i)-1:timeID(i);
%             plot(xp(dim),yp(dim),'r')
%             plot(xp(dim),yp(dim),'r*')
%         end
%         hold off
        
        cross = xinc | yinc | ydec;
        
        pdiff = ((xpt-xp1).^2 + (ypt-yp1).^2).^0.5;
        
        Vpint = Vp1 + (xpt-xp1).*((Vp2-Vp1)./(xp2-xp1));
        
        Vp2(:,cross) = Vpint(:,cross);
        
        Vavg = (Vp1+Vp2)./2;
        Vmag = (Vavg(1,:).^2 + Vavg(2,:).^2).^0.5;
        
        realtime(particleArray) = realtime(particleArray) + pdiff./Vmag;
        
        %% Stop particle
        % Stopping conditions: When particle reaches surface boundary
        maxRow = rows >= dim1;
        maxCol = cols >= dim2;
        minCol = cols <= 0;
        
        % Stop particle if any of above conditions met
        stopCon = maxRow | maxCol | minCol;
        
        % Remove such particles so that only those still within surface are
        % being calculated
        particleArray = particleArray(~stopCon);
        
        rows = rows(~stopCon);
        cols = cols(~stopCon);
        
        xCross = xCross(:,~stopCon);
        yCross = yCross(:,~stopCon);
        
        xPrevCross = xPrevCross(:,~stopCon);
        yPrevCross = yPrevCross(:,~stopCon);
        
        xp1 = xpt(:,~stopCon);
        yp1 = ypt(:,~stopCon);
        Vp1 = Vp2(:,~stopCon);
        
        % Update ID matrix for panel corner points and velocities
        ID = rows + (cols-1)*dim3;
        
        % Update time and timeID for next step
        t = t + 1;
        timeID = timeID + 1;
        timeID = timeID(~stopCon);
        
        % New corner points and velocities for next step
        x1 = ax(ID);
        x2 = bx(ID);
        x3 = cx(ID);
        x4 = dx(ID);
        
        y1 = ay(ID);
        y2 = by(ID);
        y3 = cy(ID);
        y4 = dy(ID);
        
        c1 = [x1; y1];
        c2 = [x2; y2];
        c3 = [x3; y3];
        c4 = [x4; y4];
        
        V1 = [Vx1(ID); Vy1(ID)];
        V2 = [Vx2(ID); Vy2(ID)];
        V3 = [Vx3(ID); Vy3(ID)];
        V4 = [Vx4(ID); Vy4(ID)];
        
    end
else
    
end

figure(1)
hold on
plot(xp,yp,'r')
plot(xp,yp,'r*')
hold off