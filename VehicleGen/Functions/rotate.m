function partStruct = rotate(partStruct,rotCentre,alpha)

% More accurate to rotate full configuration and calculate normals again to
% find inclination to flow. However these normals not saved in partArray as
% original normals used to find CN, CA, before transferring to CL, CD with
% angles to planes.

for i=1:length(partStruct)
    
    points = partStruct(i).Points;
    
    % Rotating in the xz plane
    points(:,:,1) = ((points(:,:,1) - rotCentre) * cos(alpha) + points(:,:,3) * sin(alpha)) + rotCentre;
    points(:,:,3) = points(:,:,3)*cos(alpha) - points(:,:,1)*sin(alpha);
    
    a = points(1:end-1,1:end-1,:);
    b = points(2:end,1:end-1,:);
    c = points(2:end,2:end,:);
    d = points(1:end-1,2:end,:);
    
    centre = (a + b + c + d)/4;
    
    %% Calculate normals of each panel
    ac = c - a;
    db = b - d;

    xNorm = db(:,:,2).*ac(:,:,3) - db(:,:,3).*ac(:,:,2);
    yNorm = db(:,:,3).*ac(:,:,1) - db(:,:,1).*ac(:,:,3);
    zNorm = db(:,:,1).*ac(:,:,2) - db(:,:,2).*ac(:,:,1);
    
    magNorm = (xNorm.^2 + yNorm.^2 + zNorm.^2).^0.5;
    
    %% CHECK
    % Zero magnitude normals (ie lines not planes) will give NaN when 
    % normalising normal components, so set to small value to avoid this
    con = magNorm == 0;
    magNorm(con) = 1e-20;
    %%
    
    norm = zeros(size(magNorm));
    
    norm(:,:,1) = xNorm;
    norm(:,:,2) = yNorm;
    norm(:,:,3) = zNorm;
    
    unitNorm = norm./magNorm;
    
    yzUnitNorm = (unitNorm(:,:,2).^2 + unitNorm(:,:,3).^2).^0.5;
    halfAngle = atan2(-unitNorm(:,:,1),yzUnitNorm)*180/pi;
    
    flow = xNorm < 0;
    
    del = round(halfAngle,10);
    
    con1 = del > 90;
    con2 = del < -90;
    del(con1) = 180 - del(con1);
    del(con2) = -180 - del(con2);
    del = del*pi/180;
    
    partStruct(i).Points = points;
    partStruct(i).norm = norm;
    partStruct(i).unitNorm = unitNorm;
    partStruct(i).centre = centre;
    partStruct(i).del = del;
    partStruct(i).flow = flow;
    
end

%% Then plot using any plotter function
% plotnorms(partArray);
% plotter(partArray);
end