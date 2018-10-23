clear
close all

% Use for incoming flow
Vinf = [100,50,0];

% Use for interpolating 4 corner plate velocities
% NOTE: Corners here defined as anti-clockwise from bottom left corner
Vc1 = [50,50];
Vc2 = [50,50];
Vc3 = [50,50];
Vc4 = [50,50];

%% Creating uniform/randomly perturbed panels
% Replace with size(points.x)
x = ones(2); % x,y points (x-1,y-1 panels)

% Creating x,y point matrices from 0 to [x y]-1
% Number of points
[dim1,dim2] = size(x);

% Number of panels
dim3 = dim1-1;
dim4 = dim2-1;

x = repmat((0:dim1-1)'*2,1,dim2);
y = repmat((0:dim2-1)*2,dim1,1);
z = zeros(dim1,dim2);

xlognorm = x(1:end-1,1:end-1);
ylognorm = y(1:end-1,1:end-1);

%% Random perturbation to make uniform grid chaotic, only inner points
% randx = rand(dim3-1,dim4-1)*0.5 - 0.25;
% randy = rand(dim3-1,dim4-1)*0.5 - 0.25;
% 
% % randx = rand(dim1-2,dim2-2) - 0.5;
% % randy = rand(dim1-2,dim2-2) - 0.5;
% 
% xr = x;
% yr = y;
% 
% xr(2:end-1,2:end-1) = x(2:end-1,2:end-1) + randx;
% yr(2:end-1,2:end-1) = y(2:end-1,2:end-1) + randy;

%% TEMP FOR ONE QUAD
randx = rand(2,2)*0.5 - 0.25;
randy = rand(2,2)*0.5 - 0.25;

xr = x/2 + randx;
yr = y/2 + randy;
%%

xyz = zeros(dim1,dim2*3);

xyz(:,1:3:end) = x;
xyz(:,2:3:end) = y;
xyz(:,3:3:end) = z;

points.Name = "test";
points.x = x;
points.y = y;
points.z = z;
points.xyz = xyz;

xyzreal = zeros(dim1,dim2*3);

xyzreal(:,1:3:end) = xr;
xyzreal(:,2:3:end) = yr;
xyzreal(:,3:3:end) = z;

physical.Name = "test";
physical.x = xr;
physical.y = yr;
physical.z = z;
physical.xyz = xyzreal;

points = normals(points);
physical = normals(physical);
plotter(points)
plotter(physical)

centre = points.centre;
centrex = centre(:,1:3:end);
centrey = centre(:,2:3:end);
norm = points.unitNorm;

nx = norm(:,1:3:end);
ny = norm(:,2:3:end);
nz = norm(:,3:3:end);

a = physical.a;
b = physical.b;
c = physical.c;
d = physical.d;

ax = a(:,1:3:end);
bx = b(:,1:3:end);
cx = c(:,1:3:end);
dx = d(:,1:3:end);

ay = a(:,2:3:end);
by = b(:,2:3:end);
cy = c(:,2:3:end);
dy = d(:,2:3:end);

% Edge mid points
chordx = (x(:,2:end) + x(:,1:end-1))/2;
chordy = (y(:,2:end) + y(:,1:end-1))/2;

spanx = (x(2:end,:) + x(1:end-1,:))/2;
spany = (y(2:end,:) + y(1:end-1,:))/2;

deltax = chordx(2:end,:) - chordx(1:end-1,:);
deltayz = spany(:,2:end) - spany(:,1:end-1);

% Initial injection points
xinject = chordx(1,:);
yinject = chordy(1,:);

%%

% Vx = [100, 95, 90, 85, 80];
% Vx = repmat(Vx',1,dim3);
% Vyz = [50, 30, 20, -10, -20];
% Vyz = repmat(Vyz,dim4,1);

%% Velocity definition

% Use for interpolating velocity from centre of panel to edge midpoints for
% uniform grids (potentially unit square?)
% [xV,yV,Vx1,Vx2,Vx3,Vx4,Vy1,Vy2,Vy3,Vy4] = centretomidvelocity(points,Vinf);

% Use for interpolating surface corner to velocities to panel corners for
% uniform grids. Both methods give basically same result

%% Simple bilinear interpolation, arbitrary rectangle turned to unit square
u = x./max(x(:)) - x(1);
v = y./max(y(:)) - y(1);

v0 = (1-u).*Vc1(1) + u.*Vc2(1);
v1 = (1-u).*Vc4(1) + u.*Vc3(1);

Vx = (1-v).*v0 + v.*v1;

v0 = (1-u).*Vc1(2) + u.*Vc2(2);
v1 = (1-u).*Vc4(2) + u.*Vc3(2);

Vy = (1-v).*v0 + v.*v1;

%% Bilinear interpolation, full form, no need to translate into unit square
mat1 = [1, xr(1), yr(1), xr(1)*yr(1);
    1, xr(1), yr(end), xr(1)*yr(end);
    1, xr(end), yr(1), xr(end)*yr(1);
    1, xr(end), yr(end), xr(end)*yr(end)];

xtemp = reshape(xr,1,[]);
ytemp = reshape(yr,1,[]);

mat2 = [ones(1,dim1*dim2); xtemp; ytemp; xtemp.*ytemp];

b = inv(mat1)'*mat2;

Vx = b(1,:)*Vc1(1) + b(2,:).*Vc4(1) + b(3,:).*Vc2(1) + b(4,:).*Vc3(1);
Vy = b(1,:)*Vc1(2) + b(2,:).*Vc4(2) + b(3,:).*Vc2(2) + b(4,:).*Vc3(2);

Vx = reshape(Vx,dim1,[]);
Vy = reshape(Vy,dim1,[]);

Vx1 = (Vx(1:end-1,1:end-1) + Vx(1:end-1,2:end))/2;
Vx2 = (Vx(2:end,1:end-1) + Vx(2:end,2:end))/2;

Vyz1 = (Vy(1:end-1,1:end-1) + Vy(2:end,1:end-1))/2;
Vyz2 = (Vy(1:end-1,2:end) + Vy(2:end,2:end))/2;

Ax = (Vx2-Vx1)/2;%./deltax;
Ayz = (Vyz2-Vyz1)/2;%./deltayz;

stopCon = false(1,dim4);

rows = ones(1,dim4);
cols = 1:dim4;

nRows = dim3*10;
dummy = nan(nRows,dim4);

xp = dummy;
yp = dummy;
xtrans = dummy;
ytrans = dummy;
Vxp = dummy;
Vyzp = dummy;

lastRow = zeros(1,dim4);

%% Initialisation

t = zeros(1,dim4);

i = 1;
% Initial matrix ID, first leading edge panels
ID = rows + (cols-1)*dim3;
timeID = i + (cols-1)*nRows;

% 1: incoming panel velocity, 2: outgoing panel velocity 
Vx1i = Vx1(ID);
Vx2i = Vx2(ID);

% Same as above but in yz plane
Vyz1i = Vyz1(ID);
Vyz2i = Vyz2(ID);

Axi = Ax(ID);
Ayzi = Ayz(ID);
% deltaxi = deltax(ID);
% deltayzi = deltayz(ID);
centrexi = centrex(ID);
centreyi = centrey(ID);

pcornx = [ax(ID); bx(ID); cx(ID); dx(ID)];
pcorny = [ay(ID); by(ID); cy(ID); dy(ID)];

% Particle injection x,y,z coordinates
xp(i,:) = xinject;
yp(i,:) = yinject;

%% Transform to real coordinates 1
% l = (xinject - xlognorm(ID));
% m = (yinject - ylognorm(ID));
% 
% % Non-dimensionalise to unit square coordinates
% %compute coefficients
% A=[1 0 0 0;1 1 0 0;1 1 1 1;1 0 1 0];
% 
% AX = A\pcornx;
% BY = A\pcorny;
% 
% xtrans(i,:) = AX(1,:) + AX(2,:).*l + AX(3,:).*m + AX(4,:).*l.*m;
% ytrans(i,:) = BY(1,:) + BY(2,:).*l + BY(3,:).*m + BY(4,:).*l.*m;

%% Transform to real coordinates 2
% eps = [-1,1,1,-1];
% eta = [-1,-1,1,1];
% 
% eps = [eps xinject - centrex(ID)];
% eta = [eta yinject - centrey(ID)];

eps = xinject - centrex(ID);
eta = yinject - centrey(ID);

N1 = 0.25*(1-eps).*(1-eta);
N2 = 0.25*(1+eps).*(1-eta);
N3 = 0.25*(1+eps).*(1+eta);
N4 = 0.25*(1-eps).*(1+eta);

xtrans(i,:) = N1.*pcornx(1,:) + N2.*pcornx(2,:) + N3.*pcornx(3,:) + N4.*pcornx(4,:);
ytrans(i,:) = N1.*pcorny(1,:) + N2.*pcorny(2,:) + N3.*pcorny(3,:) + N4.*pcorny(4,:);

dx_dxhat = 0.25*(-pcornx(1,:)+pcornx(2,:)+pcornx(3,:)-pcornx(4,:) +...
    eta.*(pcornx(1,:)-pcornx(2,:)+pcornx(3,:)-pcornx(4,:)));

dx_dyhat = 0.25*(-pcornx(1,:)-pcornx(2,:)+pcornx(3,:)+pcornx(4,:) +...
    eps.*(pcornx(1,:)-pcornx(2,:)+pcornx(3,:)-pcornx(4,:)));

dy_dxhat = 0.25*(-pcorny(1,:)+pcorny(2,:)+pcorny(3,:)-pcorny(4,:) +...
    eta.*(pcorny(1,:)-pcorny(2,:)+pcorny(3,:)-pcorny(4,:)));

dy_dyhat = 0.25*(-pcorny(1,:)-pcorny(2,:)+pcorny(3,:)+pcorny(4,:) +...
    eps.*(pcorny(1,:)-pcorny(2,:)+pcorny(3,:)-pcorny(4,:)));

D = [dx_dxhat, dx_dyhat;...
    dy_dxhat, dy_dyhat];

J = dx_dxhat.*dy_dyhat - dx_dyhat.*dy_dxhat;

invD = inv(D);

V1i = J*invD*[Vx1i; Vyz1i];
V2i = J*invD*[Vx2i; Vyz2i];

Vx1i = V1i(1);
Vyz1i = V1i(2);
Vx2i = V2i(1);
Vyz2i = V2i(2);

%% Flow direction
% Assuming positive yz direction flow to begin with
inCol = 1:dim2-1;
outCol = 2:dim2;

flip = Vyz1i < 0;
tempyz1 = Vyz2i(flip);
tempyz2 = Vyz1i(flip);
tempInCol = outCol(flip);
tempOutCol = inCol(flip);

Vyz1i(flip) = tempyz1;
Vyz2i(flip) = tempyz2;
inCol(flip) = tempInCol;
outCol(flip) = tempOutCol;

inIDx = rows + (cols-1)*dim1;
outIDx = (rows+1) + (cols-1)*dim1;

inIDy = rows + (inCol-1)*dim1;
outIDy = rows + (outCol-1)*dim1;

Vxp(timeID) = Axi.*(xp(timeID)-x(inIDx)) + Vx1i;
Vyzp(timeID) = Ayzi.*(yp(timeID)-y(inIDy)) + Vyz1i;

counter = timeID;

while any(~stopCon)
    
    xCon = Vx2i == Vxp(timeID);
    yCon = Vyz2i == Vyzp(timeID);
    
    deltatx = abs((1./Axi).*log(Vx2i./Vxp(timeID)));
    deltatyz = abs((1./Ayzi).*log(Vyz2i./Vyzp(timeID)));
    
    if any(Vx2i == Vxp(timeID))
    
        deltatx2 = (x(outIDx)-xp(timeID))./Vx1i;
        deltatx(xCon) = deltatx2(xCon);
    end
    
    if any(Vyz2i == Vyzp(timeID))
    
        deltatyz2 = (y(outIDy)-yp(timeID))./Vyz1i;
        deltatyz(yCon) = deltatyz2(yCon);
    end
    
    yspecial = (Vyz2i > 0 & Vyz1i < 0) | (Vyz2i < 0 & Vyz1i > 0);
    deltatyz(yspecial) = inf;
    
    deltat = min(deltatx,deltatyz);
    
    xpe = x(inIDx) + (1./Axi).*(Vxp(timeID).*exp(Axi.*deltat) - Vx1i);
    ype = y(inIDy) + (1./Ayzi).*(Vyzp(timeID).*exp(Ayzi.*deltat) - Vyz1i);
    
    if any(Vx2i == Vx1i)
        
        xpe2 = xp(timeID) + (deltat.*Vx1i);
        xpe(xCon) = xpe2(xCon);
    end
    
    if any(Vyz2i == Vyz1i)
    
        ype2 = yp(timeID) + (deltat.*Vyz1i);
        ype(yCon) = ype2(yCon);
    end
    
    xp(timeID+1) = xpe;
    yp(timeID+1) = ype;
    
    %% Transform to real coordinates 1 
    
%     l = xpe - xlognorm(ID);
%     m = ype - ylognorm(ID);
%     
%     AX = A\pcornx;
%     BY = A\pcorny;
%     
%     xtrans(timeID+1) = AX(1,:) + AX(2,:).*l + AX(3,:).*m + AX(4,:).*l.*m;
%     ytrans(timeID+1) = BY(1,:) + BY(2,:).*l + BY(3,:).*m + BY(4,:).*l.*m;
    
    %% Transform to real coordinates 2
    eps = xpe - centrex(ID);
    eta = ype - centrey(ID);
    
    N1 = 0.25*(1-eps).*(1-eta);
    N2 = 0.25*(1+eps).*(1-eta);
    N3 = 0.25*(1+eps).*(1+eta);
    N4 = 0.25*(1-eps).*(1+eta);
    
    xtrans(timeID+1) = N1.*pcornx(1,:) + N2.*pcornx(2,:) + N3.*pcornx(3,:) + N4.*pcornx(4,:);
    ytrans(timeID+1) = N1.*pcorny(1,:) + N2.*pcorny(2,:) + N3.*pcorny(3,:) + N4.*pcorny(4,:);
    
    %% Plotting
    figure(1)
    hold on
    for i=1:sum(~stopCon)
        dim = timeID(i):timeID(i)+1;
        plot(xp(dim),yp(dim),'r')
        plot(xp(dim),yp(dim),'r*')
    end
    hold off
    
    figure(2)
    hold on
    for i=1:sum(~stopCon)
        dim = timeID(i):timeID(i)+1;
        plot(xtrans(dim),ytrans(dim),'r')
        plot(xtrans(dim),ytrans(dim),'r*')
    end
    hold off
    
    %%
    xinc = deltatx < deltatyz;
    yinc = ~xinc & ~flip & ~yspecial;
    ydec = ~xinc & flip & ~yspecial;
    
    special = deltatx == deltatyz;
    rows(special) = rows(special) + 1;
    
    rows(xinc) = rows(xinc) + 1;
    cols(yinc) = cols(yinc) + 1;
    cols(ydec) = cols(ydec) - 1;
    
    maxRow = rows >= dim1;
    maxCol = cols >= dim2;
    minCol = cols <= 0;
    
    t(~stopCon) = t(~stopCon) + deltat;
    
    stopCon = maxRow | maxCol | minCol;
    
    rows = rows(~stopCon);
    cols = cols(~stopCon);
    
    i = i + 1;
    counter = counter + 1;
    ID = rows + (cols-1)*dim3;
    timeID = timeID + 1;
    timeID = timeID(~stopCon);
    
    Vx1i = Vx1(ID);
    Vyz1i = Vyz1(ID);
    Vx2i = Vx2(ID);
    Vyz2i = Vyz2(ID);
    
    Axi = Ax(ID);
    Ayzi = Ayz(ID);
    centrexi = centrex(ID);
    centreyi = centrey(ID);
    
    pcornx = [ax(ID); bx(ID); cx(ID); dx(ID)];
    pcorny = [ay(ID); by(ID); cy(ID); dy(ID)];
    
    inCol = cols;
    outCol = inCol+1;
    
    flip = Vyz1i < 0;
    tempyz1 = Vyz2i(flip);
    tempyz2 = Vyz1i(flip);
    tempInCol = outCol(flip);
    tempOutCol = inCol(flip);
    
    Vyz1i(flip) = tempyz1;
    Vyz2i(flip) = tempyz2;
    inCol(flip) = tempInCol;
    outCol(flip) = tempOutCol;
    
    inIDx = rows + (cols-1)*dim1;
    outIDx = (rows+1) + (cols-1)*dim1;
    
    inIDy = rows + (inCol-1)*dim1;
    outIDy = rows + (outCol-1)*dim1;
    
    Vxp(timeID) = Axi.*(xp(timeID)-x(inIDx)) + Vx1i;
    Vyzp(timeID) = Ayzi.*(yp(timeID)-y(inIDy)) + Vyz1i; 
    
end

xp = xp(1:counter,:);
yp = yp(1:counter,:);
xtrans = xtrans(1:counter,:);
ytrans = ytrans(1:counter,:);