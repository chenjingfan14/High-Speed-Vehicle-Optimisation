function partStruct = rotate(partStruct,rotCentre,alpha)

% More accurate to rotate full configuration and calculate normals again to
% find inclination to flow. However these normals not saved in partArray as
% original normals used to find CN, CA, before transferring to CL, CD with
% angles to planes.

for i=1:length(partStruct)
    
    part = partStruct(i);
    
    x = part.x - rotCentre;
    y = part.y;
    z = part.z;
    
    xrot = (x*cos(alpha) + z*sin(alpha)) + rotCentre;
    yrot = y;
    zrot = z*cos(alpha) - x*sin(alpha);
    
    [dim1,dim2] = size(x);
    
    points = zeros(dim1,dim2*3);
    
    X = 1:3:dim2*3;
    Y = X + 1;
    Z = Y + 1;
    
    points(:,X) = xrot;
    points(:,Y) = y;
    points(:,Z) = zrot;
    
    at = points(1:end-1,1:end-3);
    bt = points(2:end,1:end-3);
    ct = points(2:end,4:end);
    dt = points(1:end-1,4:end);
    
    centre = (at + bt + ct + dt)/4;
    [dim3,dim4] = size(centre);
    
    act = ct - at;
    dbt = bt - dt;
    
    X = 1:3:dim4;
    Y = X + 1;
    Z = Y + 1;
    
    acx = act(:,X);
    acy = act(:,Y);
    acz = act(:,Z);
    
    dbx = dbt(:,X);
    dby = dbt(:,Y);
    dbz = dbt(:,Z);
     
    xNorm = dby.*acz - dbz.*acy;
    yNorm = dbz.*acx - dbx.*acz;
    zNorm = dbx.*acy - dby.*acx;
    
    normnorm = (xNorm.^2 + yNorm.^2 + zNorm.^2).^0.5;
    
    % Zero magnitude normals (ie lines not planes) will give NaN when 
    % normalising normal components, so set to small value to avoid this
    con = normnorm == 0;
    normnorm(con) = 1e-20;
    
    [unitNorm,norm] = deal(zeros(dim3,dim4));
    
    xUnitNorm = xNorm./normnorm;
    yUnitNorm = yNorm./normnorm;
    zUnitNorm = zNorm./normnorm;
    
    norm(:,X) = xNorm;
    norm(:,Y) = yNorm;
    norm(:,Z) = zNorm;
    
    unitNorm(:,X) = xUnitNorm;
    unitNorm(:,Y) = yUnitNorm;
    unitNorm(:,Z) = zUnitNorm;
    
    yzUnitNorm = (yUnitNorm.^2 + zUnitNorm.^2).^0.5;
    halfAngle = atan2(-xUnitNorm,yzUnitNorm)*180/pi;
    
    flow = xNorm < 0;
    
    del = round(halfAngle,10);
    
    con1 = del > 90;
    con2 = del < -90;
    del(con1) = 180 - del(con1);
    del(con2) = -180 - del(con2);
    del = del*pi/180;
    
    partStruct(i).x = xrot;
    partStruct(i).y = yrot;
    partStruct(i).z = zrot;
    partStruct(i).xyz = points;
    partStruct(i).norm = norm;
    partStruct(i).unitNorm = unitNorm;
    partStruct(i).centre = centre;
    partStruct(i).del = del;
    partStruct(i).flow = flow;
    
end

%% Then plot using any plotter function
% plotnorms(partArray);
% plotter(partArray);
