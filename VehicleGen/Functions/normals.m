function partStruct = normals(partStruct)
%% Additional panel specific properties
% Calculate centre, outward facing normals, area

for ii=1:length(partStruct)
    
    % Initialise and calculate centre points
    x = partStruct(ii).x;
    y = partStruct(ii).y;
    z = partStruct(ii).z;
    points = partStruct(ii).xyz;
    
    [dim1,dim2] = size(x);
    
    dim1 = dim1 - 1;
    dim2 = dim2 - 1;
    dim3 = dim2*3;
    
    at = points(1:end-1,1:end-3);
    bt = points(2:end,1:end-3);
    ct = points(2:end,4:end);
    dt = points(1:end-1,4:end);
    
    centre = (at + bt + ct + dt)/4;
    
    X = 1:3:dim3;
    Y = X + 1;
    Z = Y + 1;
    
    yCentre = centre(:,Y);
    zCentre = centre(:,Z);
    
    radialLocation = atan2d(yCentre - y(1),zCentre - z(1));
    
    %% Calculate normals of each panel
    act = ct - at;
    dbt = bt - dt;
    
    acx = act(:,X);
    acy = act(:,Y);
    acz = act(:,Z);
    
    dbx = dbt(:,X);
    dby = dbt(:,Y);
    dbz = dbt(:,Z);
    
    xNorm = dby.*acz - dbz.*acy;
    yNorm = dbz.*acx - dbx.*acz;
    zNorm = dbx.*acy - dby.*acx;
    
    % Magnitude of normal
    magNorm = (xNorm.^2 + yNorm.^2 + zNorm.^2).^0.5;

    [norm,unitNorm] = deal(zeros(dim1,dim3));
    
    norm(:,X) = xNorm;
    norm(:,Y) = yNorm;
    norm(:,Z) = zNorm;
    
    nx = xNorm./magNorm;
    ny = yNorm./magNorm;
    nz = zNorm./magNorm;

    unitNorm(:,X) = nx;
    unitNorm(:,Y) = ny;
    unitNorm(:,Z) = nz;
    
    % Some normal magnitudes are zero which leads to NaN. Replace such
    % values with zero
    con = isnan(unitNorm);
    unitNorm(con) = 0;
    
    nyz = (ny.^2 + nz.^2).^0.5;
    halfAngle = atan2(-nx,nyz)*180/pi;
    
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
    deltax = (x(2:end,:)-x(1:end-1,:)).^2 + (y(2:end,:)-y(1:end-1,:)).^2 + ...
        (z(2:end,:)-z(1:end-1,:)).^2;
    
    deltay = (x(:,2:end)-x(:,1:end-1)).^2 + (y(:,2:end)-y(:,1:end-1)).^2 + ...
        (z(:,2:end)-z(:,1:end-1)).^2;
    
    partStruct(ii).a = at;
    partStruct(ii).b = bt;
    partStruct(ii).c = ct;
    partStruct(ii).d = dt;
    partStruct(ii).centre = centre;
    partStruct(ii).radialLocation = radialLocation;
    partStruct(ii).norm = norm;
    partStruct(ii).unitNorm = unitNorm;
    partStruct(ii).deltax = deltax;
    partStruct(ii).deltay = deltay;
    partStruct(ii).area = panelarea(x,y,z);
    partStruct(ii).del = del;
    partStruct(ii).flow = flow;
    
end

end