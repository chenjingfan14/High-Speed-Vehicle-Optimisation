function partStruct = normals(points,viscous,quad)
%% Additional panel specific properties
% Calculate centre, outward facing normals, area

% Initialise and calculate centre points

[n,m] = size(points(:,:,1));
mn = n * m;

a = (1:mn - n)';

delete = n:n:mn - n;
a(delete,:) = [];

quadID = [a, a+1, n+a, n+a+1];

[dim,~,~] = size(a);

triID = zeros(dim*2,3);
triID(1:2:end,:) = quadID(:,1:3);
triID(2:2:end,:) = quadID(:,2:4);

x = points(:,:,1);
y = points(:,:,2);
z = points(:,:,3);

m = m - 1;
n = n - 1;

if quad
    
    a = points(1:end-1,1:end-1,:);
    b = points(2:end,1:end-1,:);
    c = points(2:end,2:end,:);
    d = points(1:end-1,2:end,:);
    
    centre = (a + b + c + d)/4;
        
    vec1 = c - a;
    vec2 = b - d;
    
    trans = size(centre);
    
    norm = crossmat(vec2, vec1);
    
    ID = reshape(1:m*n,n,m);
    
    partStruct.QuadID = quadID;
else
    xyz(:,:,1) = x(triID);
    xyz(:,:,2) = y(triID);
    xyz(:,:,3) = z(triID);
    
    trans = [2*n m 3];
    
    centre = mean(xyz,2);
    centre = reshape(centre,trans);
    
    vec1 = xyz(:,2,:) - xyz(:,1,:);
    vec2 = xyz(:,3,:) - xyz(:,1,:);
    
    vec1 = reshape(vec1,trans);
    vec2 = reshape(vec2,trans);

    norm = zeros(trans);
    
    norm(1:2:end,:,:) = crossmat(vec1(1:2:end,:,:), vec2(1:2:end,:,:));
    norm(2:2:end,:,:) = crossmat(vec2(2:2:end,:,:), vec1(2:2:end,:,:));
    
    ID = reshape(1:m*n*2,n*2,m);
end
% 
radialLocation = atan2d(centre(:,:,2) - points(1,1,2), centre(:,:,3) - points(1,1,3));

%% Calculate normals of each panel
% Normal vector and its magnitude
magNorm = magmat(norm);

unitNorm = norm./magNorm;

% Some normal magnitudes are zero which leads to NaN. Replace such
% values with zero
unitNorm(isnan(unitNorm)) = 0;
% 
% nx = unitNorm(:,:,1);
% nyz = (unitNorm(:,:,2).^2 + unitNorm(:,:,3).^2).^0.5;
% halfAngle = atan2(-nx, nyz) * 180/pi;
% 
% flow = nx < 0;
% 
% del = round(halfAngle,10);
% 
% con1 = del > 90;
% con2 = del < -90;
% del(con1) = 180 - del(con1);
% del(con2) = -180 - del(con2);
% del = del*pi/180;

%% Calculate average deltas across each panel

% Might not be needed if not using Pollock's method for viscous
% particle tracking
% xDiff = diff(points,1,1);
% yDiff = diff(points,1,2);
% 
% deltax = (xDiff(:,:,1).^2 + xDiff(:,:,2).^2 + xDiff(:,:,3).^2).^0.5;
% deltay = (yDiff(:,:,1).^2 + yDiff(:,:,2).^2 + yDiff(:,:,3).^2).^0.5;

%% Calculate panel areas
% Areas based on xyz vectors so any triangle can be input

ons = [1 1 1];

x = points(:,:,1);
y = points(:,:,2);
z = points(:,:,3);

for i = dim:-1:1
    
    I1 = i*2-1;
    I2 = I1 + 1;
    
    t1ID = triID(I1,:);
    t2ID = triID(I2,:);
    
    x2 = x(t2ID);
    y2 = y(t2ID);
    z2 = z(t2ID);
    
    A2 = 0.5*sqrt(det([x2; y2; ons])^2 + det([y2; z2; ons])^2 + det([z2; x2; ons])^2);
    
    xyzTri(I2,:,1) = x2;
    xyzTri(I2,:,2) = y2;
    xyzTri(I2,:,3) = z2;
    
    x1 = x(t1ID);
    y1 = y(t1ID);
    z1 = z(t1ID);
    
    A1 = 0.5*sqrt(det([x1; y1; ons])^2 + det([y1; z1; ons])^2 + det([z1; x1; ons])^2);
    
    xyzTri(I1,:,1) = x1;
    xyzTri(I1,:,2) = y1;
    xyzTri(I1,:,3) = z1;
    
    if ~any([A1,A2])
        
        error('Area fix required')
    end
    
    if quad
        
        At1(i,:) = A1;
        At2(i,:) = A2;
        
        area(i,:) = A1 + A2;
    else
        area(I2,:) = A2;
        area(I1,:) = A1;
    end
end

area = reshape(area,trans([1 2]));

if viscous && quad
    
    triArea = zeros(dim*2,1);
    triArea(1:2:end) = At1(:);
    triArea(2:2:end) = At2(:);
    
    triCentre = mean(xyzTri,2);

    u = xyzTri(:,2,:) - xyzTri(:,1,:);
    v = xyzTri(:,3,:) - xyzTri(:,1,:);

    triNorm = zeros(size(triCentre));
    triNorm(1:2:end,:,:) = crossmat(u(1:2:end,:,:), v(1:2:end,:,:));
    triNorm(2:2:end,:,:) = crossmat(v(2:2:end,:,:), u(2:2:end,:,:));

    magTriNorm = magmat(triNorm);

    unitTriNorm = triNorm./magTriNorm;
    
    unitTriNorm(isnan(unitTriNorm)) = 0;
    
    triangle.TriID = triID;
    triangle.Points = xyzTri;
    triangle.centre = triCentre;
    triangle.norm = triNorm;
    triangle.unitNorm = unitTriNorm;
    triangle.Area = triArea;    
    partStruct.Triangle = triangle;
end

yzPos = norm(:,:,3) > 0;

partStruct.ID = ID;
partStruct.TriID = triID;
partStruct.Points = points;
partStruct.centre = centre;
partStruct.norm = norm;
partStruct.unitNorm = unitNorm;
partStruct.area = area;
partStruct.zPos = yzPos;
% partStruct.flow = flow;
partStruct.radialLocation = radialLocation;
% partStruct.a = a;
% partStruct.b = b;
% partStruct.c = c;
% partStruct.d = d;
% partStruct.deltax = deltax;
% partStruct.deltay = deltay;

% plotter(triangle,"triangle","centre","unit normal")
% plotter(partStruct,"triangle","centre","unit normal")
end