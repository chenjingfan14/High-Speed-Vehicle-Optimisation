function [xV,yV,Vx1,Vx2,Vx3,Vx4,Vy1,Vy2,Vy3,Vy4] = centretomidvelocity(points,Vinf)

x = points.x;
y = points.y;
z = points.z;

[dim1,dim2] = size(x);

norm = points.unitNorm;
nx = norm(:,1:3:end);
ny = norm(:,2:3:end);
nz = norm(:,3:3:end);

dim3 = dim1-1;
dim4 = dim2-1;

centrex = points.centre(:,1:3:end);
centrey = points.centre(:,2:3:end);
centrez = points.centre(:,3:3:end);

%% Mid/Injection points

x13 = (x(:,2:end) + x(:,1:end-1))/2;
y13 = (y(:,2:end) + y(:,1:end-1))/2;

x24 = (x(2:end,:) + x(1:end-1,:))/2;
y24 = (y(2:end,:) + y(1:end-1,:))/2;

xV = x13(1:end-1,:);
yV = y24(:,1:end-1);

%% Surface velocity at panel centre

Vx = Vinf(1);
Vy = Vinf(2);
Vz = Vinf(3);

Tx = ny.*Vz - nz.*Vy;
Ty = nz.*Vx - nx.*Vz;
Tz = nx.*Vy - ny.*Vx;

Sx = Ty.*nz - Tz.*ny;
Sy = Tz.*nx - Tx.*nz;
Sz = Tx.*ny - Ty.*nx;

%% Mid/Corner point velocity interpolation
% Do u and v correlate to corner or mid??

Sxi = zeros(dim1,dim4);
Syi = zeros(dim1,dim4);

i = 1;
for count=1:dim1

    xi = Sx(i,:)+(x13(count,:)-centrex(i,:)).*((Sx(i+1,:)-Sx(i,:))./(centrex(i+1,:)-centrex(i,:)));
    yi = Sy(i,:)+(y13(count,:)-centrey(i,:)).*((Sy(i+1,:)-Sy(i,:))./(centrey(i+1,:)-centrey(i,:)));
%     zi(i,:) = Sz(j,:)+(z13(i,:)-centrex(j,:)).*((Sz(j+1,:)-Sz(j,:))./(centrez(j+1,:)-centrez(j,:)));
    
    % Can only = NaN if centre(i,...) = centre(i+1,...) = x/y/z13(count,...)
    % Therefore will equal either centre velocity
    reset = isnan(xi);
    xi(reset) = Sx(i,reset);

    reset = isnan(yi);
    yi(reset) = Sy(i,reset);
    
    Sxi(count,:) = xi;
    Syi(count,:) = yi;
    
    if ~any(count == [1,dim1-1])
        i = i+1;
    end
    
end

Sxj = zeros(dim3,dim2);
Syj = zeros(dim3,dim2);

j = 1;
for count=1:dim2

    xj = Sx(:,j)+(x24(:,count)-centrex(:,j)).*((Sx(:,j+1)-Sx(:,j))./(centrex(:,j+1)-centrex(:,j)));
    yj = Sy(:,j)+(y24(:,count)-centrey(:,j)).*((Sy(:,j+1)-Sy(:,j))./(centrey(:,j+1)-centrey(:,j)));
%     zi(i,:) = Sz(j,:)+(z13(i,:)-centrex(j,:)).*((Sz(j+1,:)-Sz(j,:))./(centrez(j+1,:)-centrez(j,:)));
    
    reset = isnan(xj);
    xj(reset) = Sx(reset,j);

    reset = isnan(yj);
    yj(reset) = Sy(reset,j);
    
    Sxj(:,count) = xj;
    Syj(:,count) = yj;
    
    if ~any(count == [1,dim2-1])
        j = j+1;
    end
    
end

% Vx1 = Sxi(1:end-1,:);
% Vx2 = Sxj(:,2:end);
% Vx3 = Sxi(2:end,:);
% Vx4 = Sxj(:,1:end-1);
% 
% Vy1 = Syi(1:end-1,:);
% Vy2 = Syj(:,2:end);
% Vy3 = Syi(2:end,:);
% Vy4 = Syj(:,1:end-1);

Vx1 = Sxi(1:end-1,:);
Vx2 = Sxj(:,1:end-1);
Vx3 = Sxi(2:end,:);
Vx4 = Sxj(:,2:end);

Vy1 = Syi(1:end-1,:);
Vy2 = Syj(:,1:end-1);
Vy3 = Syi(2:end,:);
Vy4 = Syj(:,2:end);