function partArray = rotate(partArray,rotCentre,alpha)

% More accurate to rotate full configuration and calculate normals again to
% find inclination to flow. However these normals not saved in partArray as
% original normals used to find CN, CA, before transferring to CL, CD with
% angles to planes.

dim = numel(partArray);

for i = dim:-1:1
    
    part = partArray(i);
    
    x = part.x - rotCentre;
    y = part.y;
    z = part.z;
    
    cx = part.centre(:,1:3:end) - rotCentre;
    cy = part.centre(:,2:3:end);
    cz = part.centre(:,3:3:end);
    
    xNorm = part.norm(:,1:3:end);
    yNorm = part.norm(:,2:3:end);
    zNorm = part.norm(:,3:3:end);
    
    xUnitNorm = part.unitNorm(:,1:3:end);
    yUnitNorm = part.unitNorm(:,2:3:end);
    zUnitNorm = part.unitNorm(:,3:3:end);
    
    xRot = (x*cos(alpha) + z*sin(alpha)) + rotCentre;
    yRot = y;
    zRot = z*cos(alpha) - x*sin(alpha);
    
    [dim1,dim2] = size(x);
    dim3 = dim1 - 1;
    dim4 = dim2 - 1;
    
    points = zeros(dim1,dim2*3);
    
    X = 1:3:dim2*3;
    Y = X + 1;
    Z = Y + 1;
    
    points(:,X) = xRot;
    points(:,Y) = y;
    points(:,Z) = zRot;
    
    a = points(1:end-1,1:end-3);
    b = points(2:end,1:end-3);
    c = points(2:end,4:end);
    d = points(1:end-1,4:end);
    
    cxRot = (cx*cos(alpha) + cz*sin(alpha)) + rotCentre;
    cyRot = cy;
    czRot = cz*cos(alpha) - cx*sin(alpha);
    
    X = 1:3:dim4*3;
    Y = X + 1;
    Z = Y + 1;
    
    [centre,normRot,unitNormRot] = deal(zeros(dim3,dim4*3));
    
    centre(:,X) = cxRot;
    centre(:,Y) = cyRot;
    centre(:,Z) = czRot;
    
    xUnitNormRot = (xUnitNorm*cos(alpha) + zUnitNorm*sin(alpha));
    yUnitNormRot = yUnitNorm;
    zUnitNormRot = zUnitNorm*cos(alpha) - xUnitNorm*sin(alpha);
    
    xNormRot = (xNorm*cos(alpha) + zNorm*sin(alpha));
    yNormRot = yNorm;
    zNormRot = zNorm*cos(alpha) - xNorm*sin(alpha);
    
    normRot(:,X) = xNormRot;
    normRot(:,Y) = yNormRot;
    normRot(:,Z) = zNormRot;
    
    unitNormRot(:,X) = xUnitNormRot;
    unitNormRot(:,Y) = yUnitNormRot;
    unitNormRot(:,Z) = zUnitNormRot;
    
    yzUnitNormRot = (yUnitNormRot.^2 + zUnitNormRot.^2).^0.5;
    halfAngle = atan2(-xUnitNormRot,yzUnitNormRot)*180/pi;
    
    flow = xUnitNorm < 0;
    
    del = round(halfAngle,10);
    
    con1 = del > 90;
    con2 = del < -90;
    del(con1) = 180 - del(con1);
    del(con2) = -180 - del(con2);
    del = del*pi/180;
    
    partArray(i).x = xRot;
    partArray(i).y = yRot;
    partArray(i).z = zRot;
    partArray(i).xyz = points;
    partArray(i).a = a;
    partArray(i).b = b;
    partArray(i).c = c;
    partArray(i).d = d;
    partArray(i).unitNorm = unitNormRot;
    partArray(i).norm = normRot;
    partArray(i).centre = centre;
    partArray(i).del = del;
    partArray(i).flow = flow;
    
end

end

%% Then plot using any plotter function
% plotnorms(partArray);
% plotter(partArray);
