function [aerofoil,aftBody,success,reason] = bodyfoil(aftBody,aerofoil,offset)

%% Get body normals for line-plane intersection (wing-body intersection)
% Bring in point matrix and x,y,z coordinates of body
bodyPoints = aftBody.Points.xyz;
xBody = aftBody.Points.x(:,1);
yBody = aftBody.Points.y(1,:);
zBody = aftBody.Points.z(1,:);

% Creating vectors for cross product to get body normals
[~,ydim] = size(yBody);
[vec1,vec2] = deal(zeros(ydim-1,3));

for j=1:ydim-1
    col = (j*3)-2:j*3;
    vec1(j,:) = bodyPoints(2:end,col+3)-bodyPoints(1:end-1,col);
    vec2(j,:) = bodyPoints(2:end,col)-bodyPoints(1:end-1,col+3);
end

% Find normals
norm = cross(vec1,vec2);

%% Bridge points
% Total number of aerofoils
numAerofoils = length(aerofoil);

[~,bDim] = size(zBody);

% Loop to find all wing-body bridges
for i=1:numAerofoils
    
    wingtail = aerofoil(i);
    
    bool = wingtail.Boolean;
    foil = wingtail.Points.xyz;
    
    % Apply foil x-offset wrt aft-body
    aftOffset = offset(1) + xBody(1);
    foil(:,1:3:end) = foil(:,1:3:end) - foil(1,1) + aftOffset;
    foil(:,3:3:end) = foil(:,3:3:end) - foil(1,3) + offset(2);
    
    zFoil = foil(:,3:3:end);
    
    % Body points between min and max z
    ind = 1:bDim;
    cond = zBody >= min(zFoil(:)) & zBody <= max(zFoil(:));
    
    if any(cond)
        
        ind(~cond) = [];

        % Upper and lower boundaries of radial body panels to include in plane
        % intersection calculation ie. potential panels that bridge point could
        % lie on. Wing not allowed to interfere with first or last radial
        % panels hence max(...,2), min(...,bDim-2) (Should be bDim-1?)
        % Pull in additional points just outside min (-1) and max (+1) z
        ub = max(min(ind)-1,2);
        lb = min(max(ind)+1,bDim-2);
        
    else
        % Upper/Lower boundaries may be entirely within two radial body
        % panels, thus no body points between min and max for cond. So
        % use all radial panels (bar first and last)
        ub = 2;
        lb = bDim-2;
        
    end
    
    j = ub:lb;
    col = (ub*3)-2:lb*3;
    
    interiorPoints = [foil(:,1:3);foil(:,end-2:end)];
    exteriorPoints = [foil(:,4:6);foil(:,end-5:end-3)];
    
    [fRow,fCol] = size(interiorPoints);
    bridge = zeros(fRow,fCol);
    
    % Calculate bridge points ie. point on body where aerofoil joins
    for ii=1:fRow
        interiorPoint = interiorPoints(ii,:);
        exteriorPoint = exteriorPoints(ii,:);
        
        [inter,reason] = planeintersection(norm(j,:),bodyPoints(end,col)',interiorPoint,exteriorPoint,xBody,yBody',zBody',j);
        if isempty(inter)
            
            % Check what failed config looks like
            aerofoil.Points.x = foil(:,1:3:end);
            aerofoil.Points.z = foil(:,3:3:end);
            plotter([aftBody.Points,aerofoil.Points]);
            
            success = false;
            return
        end
        bridge(ii,:) = inter;
    end
    
    if bool
        lowerBridge = bridge;
        upperBridge = zeros(size(bridge));
    else
        fRow = fRow/2;
        upperBridge = bridge(1:fRow,:);
        lowerBridge = bridge(fRow+1:end,:);
    end
    
    % Finds where other aerofoils interfere with said panels and futher discretises
    % only these panels
    xUpperBridge = upperBridge(:,1);
    xLowerBridge = lowerBridge(:,1);
    
    yUpperBridge = upperBridge(:,2);
    yLowerBridge = lowerBridge(:,2);
    
    zUpperBridge = upperBridge(:,3);
    zLowerBridge = lowerBridge(:,3);
    
    LE = bridge(1,:); % Leading edge z point
    
    zAbsDiff = abs(zBody - LE(3));
    [~,upCut] = min(zAbsDiff);
    loCut = bDim-upCut;
    
    xUpperBridgeMat = repmat(xUpperBridge,1,upCut);
    xLowerBridgeMat = repmat(xLowerBridge,1,loCut+1);
    
    yUpperBridgeMat = repmat(yUpperBridge,1,upCut);
    yLowerBridgeMat = repmat(yLowerBridge,1,loCut+1);
    
    zUpperBridgeMat = repmat(zUpperBridge,1,upCut);
    zLowerBridgeMat = repmat(zLowerBridge,1,loCut+1);
    
    idUpper = 1:upCut;
    idLower = (0:loCut) + upCut;
    
    yBodyUpper = repmat(yBody(idUpper),fRow,1);
    yBodyLower = repmat(yBody(idLower),fRow,1);
    
    zBodyUpper = repmat(zBody(idUpper),fRow,1);
    zBodyLower = repmat(zBody(idLower),fRow,1);
    
    conUpper = zBodyUpper < zUpperBridgeMat;
    conLower = zBodyLower > zLowerBridgeMat;
    
    % Set closest bodypoint lines equal to wing bridges
    conUpper(:,end) = true;
    conLower(:,1) = true;
    
    xBodyUpper = xUpperBridgeMat;
    xBodyLower = xLowerBridgeMat;
    
    yBodyUpper(conUpper) = yUpperBridgeMat(conUpper);
    yBodyLower(conLower) = yLowerBridgeMat(conLower);
    
    zBodyUpper(conUpper) = zUpperBridgeMat(conUpper);
    zBodyLower(conLower) = zLowerBridgeMat(conLower);
    
    %% Adding stupid extra beginning end points
    
    % TODO: Put in loop
    % z(end) = z(end-1) to solve collapsing trailing edges
    
    xu = repmat(xBody,1,upCut);
    xl = repmat(xBody,1,loCut+1);
    
%     aftBody.Points(1).Name = "aftbody";
    aftBody.Points(1).x = [xu(1,:); xBodyUpper; xu(end,:)];
    aftBody.Points(1).y = [yBodyUpper(1,:); yBodyUpper; yBodyUpper(end,:)];
    aftBody.Points(1).z = [zBodyUpper(1,:); zBodyUpper; zBodyUpper(end,:)];
    
%     aftBody.Points(2).Name = "aftbody";
    aftBody.Points(2).x = [xl(1,:); xBodyLower; xl(end,:)];
    aftBody.Points(2).y = [yBodyLower(1,:); yBodyLower; yBodyLower(end,:)];
    aftBody.Points(2).z = [zBodyLower(1,:); zBodyLower; zBodyLower(end,:)];
    
    % Wing first and last point sets now replaced with upper and lower
    % bridge points found above
    wingtail.Points.xyz = [upperBridge,foil(:,4:end-3),lowerBridge];
    wingtail.Points = pointstoxyz(wingtail.Points);
    
    wetChord1 = diff(xUpperBridge([1 end]));
    wetChord2 = wingtail.Chord(2);
    
    y1 = reshape(wingtail.Points.y(:,[1 end]),[],1);
    y2 = reshape(wingtail.Points.y(:,[2 end-1]),[],1);
    wetSpan = diff(mean([y1 y2],1));
    
    wetArea = ((wetChord2 + wetChord1)/2)*wetSpan;
    Taper = wetChord2/wetChord1;
    if Taper == inf
        Taper = 0;
    end
    wetMAC = (2/3)*wetChord1*((1+Taper+(Taper^2))/(1+Taper));
    
    wingtail.WetChord = [diff(xUpperBridge([1 end])), wingtail.Chord(2:end)];
    wingtail.WetArea = [wetArea, wingtail.Area(2:end)];
    wingtail.WetSpan = [wetSpan, wingtail.Span(2:end)];
    wingtail.WetMAC = [wetMAC, wingtail.MAC(2:end)];
    
    aerofoil(i) = wingtail;
end

aftBody.Points = xyztopoints(aftBody.Points);
success = true;