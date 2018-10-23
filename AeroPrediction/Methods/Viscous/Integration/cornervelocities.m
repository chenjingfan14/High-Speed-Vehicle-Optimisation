function points = cornervelocities(points,rot,flow)

% clear all
% close all
% 
% flow = flowparameters();
% points = plategen();

[~,parts] = size(points);

Vinf = flow.U;

for i = 1:parts
    
    part = points(i);
    
    x = part.x;
    y = part.y;
    z = part.z;
    xyz = part.xyz;
    
    [dim1,dim2] = size(xyz);
    
    dim3 = dim1-1;
    dim4 = dim2-3;
    
    X = 1:3:dim2;
    Y = X + 1;
    Z = Y + 1;
    
    Vc = zeros(dim1,dim2);
    
    r1 = 1:dim1;
    c1 = 1:dim2;
    
    r2 = 1:dim3;
    c2 = 1:dim4;
    innerRows = (2:dim3)';
    innerCols = (4:dim4);
    
    ID = (1:dim1)' + ((1:dim2)-1)*dim1;
    
    a = [1,r2]' + (c1-1)*dim1;
    b = r1' + ([1:3,c2]-1)*dim1;
    c = [r2,dim3]'+1 + (c1-1)*dim1;
    d = r1' + ([c2,dim4-2:dim4]+2)*dim1;
    
    pa = xyz(a) - xyz;
    pb = xyz(b) - xyz;
    pc = xyz(c) - xyz;
    pd = xyz(d) - xyz;
    
    pax = pa(:,X);
    pay = pa(:,Y);
    paz = pa(:,Z);
    
    pbx = pb(:,X);
    pby = pb(:,Y);
    pbz = pb(:,Z);
    
    pcx = pc(:,X);
    pcy = pc(:,Y);
    pcz = pc(:,Z);
    
    pdx = pd(:,X);
    pdy = pd(:,Y);
    pdz = pd(:,Z);
    
    xNorm1 = pay.*pbz - paz.*pby;
    yNorm1 = paz.*pbx - pax.*pbz;
    zNorm1 = pax.*pby - pay.*pbx;
    
    xNorm2 = pby.*pcz - pbz.*pcy;
    yNorm2 = pbz.*pcx - pbx.*pcz;
    zNorm2 = pbx.*pcy - pby.*pcx;
    
    xNorm3 = pcy.*pdz - pcz.*pdy;
    yNorm3 = pcz.*pdx - pcx.*pdz;
    zNorm3 = pcx.*pdy - pcy.*pdx;
    
    xNorm4 = pdy.*paz - pdz.*pay;
    yNorm4 = pdz.*pax - pdx.*paz;
    zNorm4 = pdx.*pay - pdy.*pax;
    
    magNorm1 = (xNorm1.^2 + yNorm1.^2 + zNorm1.^2).^0.5;
    magNorm2 = (xNorm2.^2 + yNorm2.^2 + zNorm2.^2).^0.5;
    magNorm3 = (xNorm3.^2 + yNorm3.^2 + zNorm3.^2).^0.5;
    magNorm4 = (xNorm4.^2 + yNorm4.^2 + zNorm4.^2).^0.5;
    
    con1 = magNorm1 > 0;
    con2 = magNorm2 > 0;
    con3 = magNorm3 > 0;
    con4 = magNorm4 > 0;
    
    numNorms = con1 + con2 + con3 + con4;
    
    xNorm = (xNorm1 + xNorm2 + xNorm3 + xNorm4)./numNorms;
    yNorm = (yNorm1 + yNorm2 + yNorm3 + yNorm4)./numNorms;
    zNorm = (zNorm1 + zNorm2 + zNorm3 + zNorm4)./numNorms;
    
    magNorm = (xNorm.^2 + yNorm.^2 + zNorm.^2).^0.5;
    
    con = isnan(magNorm);
    
    % Check signs for positive z etc
    if any(con)
        xNorm(con) = -cos(rot);
        yNorm(con) = 0;
        zNorm(con) = -sin(rot);
        
        magNorm(con) = (xNorm(con).^2 + yNorm(con).^2 + zNorm(con).^2).^0.5;
    end
    
    nx = xNorm./magNorm;
    ny = yNorm./magNorm;
    nz = zNorm./magNorm;
    
%     nx = part.unitNorm(:,1:3:end);
%     ny = part.unitNorm(:,2:3:end);
%     nz = part.unitNorm(:,3:3:end);

    xZero = isnan(nx);
    yZero = isnan(ny);
    zZero = isnan(nz);
    
    nx(xZero) = 0;
    ny(yZero) = 0;
    nz(zZero) = 0;
    
%     figure(1)
%     hold on
%     axis equal
%     plot3(x+xNorm,y+yNorm,z+zNorm,'r*')
%     hold off
    
    %% Centre point velocities
    Vx = Vinf(1);
    Vy = Vinf(2);
    Vz = Vinf(3);
    
    Tx = ny.*Vz - nz.*Vy;
    Ty = nz.*Vx - nx.*Vz;
    Tz = nx.*Vy - ny.*Vx;
    
    Sx = Ty.*nz - Tz.*ny;
    Sy = Tz.*nx - Tx.*nz;
    Sz = Tx.*ny - Ty.*nx;
    
    Vc(:,X) = Sx;
    Vc(:,Y) = Sy;
    Vc(:,Z) = Sz;
    
    points(i).Vc = Vc;

end