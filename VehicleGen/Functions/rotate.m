function partStruct = rotate(partStruct,rotCentre,alpha,quad)

% More accurate to rotate full configuration and calculate normals again to
% find inclination to flow. However these normals not saved in partArray as
% original normals used to find CN, CA, before transferring to CL, CD with
% angles to planes.

for i=1:length(partStruct)
    
    points = partStruct(i).Points;
    centre = partStruct(i).centre;
    
    [n,m,~] = size(points);
    
    % Rotating in the xz plane
    points(:,:,1) = ((points(:,:,1) - rotCentre) * cos(alpha) + points(:,:,3) * sin(alpha)) + rotCentre;
    points(:,:,3) = points(:,:,3)*cos(alpha) - points(:,:,1)*sin(alpha);
    
    centre(:,:,1) = ((centre(:,:,1) - rotCentre) * cos(alpha) + centre(:,:,3) * sin(alpha)) + rotCentre;
    centre(:,:,3) = centre(:,:,3)*cos(alpha) - centre(:,:,1)*sin(alpha);
    
    if quad
        
        a = points(1:end-1,1:end-1,:);
        b = points(2:end,1:end-1,:);
        c = points(2:end,2:end,:);
        d = points(1:end-1,2:end,:);
        
        vec1 = c - a;
        vec2 = b - d;
        
        norm = crossmat(vec2, vec1);
    else
        triID = partStruct(i).TriID;
        
        x = points(:,:,1);
        y = points(:,:,2);
        z = points(:,:,3);
        
        xyz = zeros((n-1)*(m-1)*2,3,3);
        
        xyz(:,:,1) = x(triID);
        xyz(:,:,2) = y(triID);
        xyz(:,:,3) = z(triID);
        
        trans = [2*(n-1) m-1 3];
        
        vec1 = xyz(:,2,:) - xyz(:,1,:);
        vec2 = xyz(:,3,:) - xyz(:,1,:);
        
        vec1 = reshape(vec1,trans);
        vec2 = reshape(vec2,trans);
        
        norm = zeros(trans);
        
        norm(1:2:end,:,:) = crossmat(vec1(1:2:end,:,:), vec2(1:2:end,:,:));
        norm(2:2:end,:,:) = crossmat(vec2(2:2:end,:,:), vec1(2:2:end,:,:));
    end
    
    magNorm = magmat(norm);
    
    unitNorm = norm./magNorm;
    unitNorm(isnan(unitNorm)) = 0;
    
%     yzUnitNorm = (unitNorm(:,:,2).^2 + unitNorm(:,:,3).^2).^0.5;
%     halfAngle = atan2(-unitNorm(:,:,1),yzUnitNorm)*180/pi;
    
    flow = norm(:,:,1) < 0;
    
%     del = round(halfAngle,10);
%     
%     con1 = del > 90;
%     con2 = del < -90;
%     del(con1) = 180 - del(con1);
%     del(con2) = -180 - del(con2);
%     del = del*pi/180;
    
    partStruct(i).Points = points;
    partStruct(i).norm = norm;
    partStruct(i).unitNorm = unitNorm;
    partStruct(i).centre = centre;
%     partStruct(i).del = del;
    partStruct(i).flow = flow;
end

%% Then plot using any plotter function
% plotnorms(partArray);
% plotter(partArray);
end