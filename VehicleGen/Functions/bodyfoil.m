function [aerofoil,aftBody,success,reason] = bodyfoil(aftBody,aerofoil,offset)

%% Get body normals for line-plane intersection (wing-body intersection)
% Bring in point matrix and x,y,z coordinates of body
bodyPoints = aftBody.Points;

% Chordwise x body coords
xBody = bodyPoints(:,1,1);
% Radial y,z body coords
yzBody = squeeze(bodyPoints(1,:,[2 3]));

% Creating vectors for cross product to get body normals
[row,col,~] = size(bodyPoints);
rows = 1:row - 1;
cols = 1:col - 1;

a = bodyPoints(rows,cols,:);
b = bodyPoints(rows + 1,cols,:);
c = bodyPoints(rows + 1,cols + 1,:);
d = bodyPoints(rows,cols + 1,:);

ac = c - a;
db = b - d;

% Find normals. Indexing will only work for single panel, constant shape
% aftbody
norm(:,1) = db(:,:,2).*ac(:,:,3) - db(:,:,3).*ac(:,:,2);
norm(:,2) = db(:,:,3).*ac(:,:,1) - db(:,:,1).*ac(:,:,3);
norm(:,3) = db(:,:,1).*ac(:,:,2) - db(:,:,2).*ac(:,:,1);

norm = -norm;

%% Bridge points
% Total number of aerofoils
numAerofoils = length(aerofoil);

% Loop to find all wing-body bridges
for i=1:numAerofoils
    
    wingtail = aerofoil(i);
    
    bool = wingtail.Boolean;
    foil = wingtail.Points;
    
    % Remove any offsets from aerofoil
    foil = foil - foil(1,1,:);
    
    % Apply foil x-offset wrt aft-body
    aftOffset = offset(1) + xBody(1);
    foil(:,:,1) = foil(:,:,1) + aftOffset;
    foil(:,:,3) = foil(:,:,3) + offset(2);
    
    zFoil = foil(:,:,3);
    
    % Body points between min and max z
    ind = 1:col;
    cond = yzBody(:,2) >= min(zFoil(:)) & yzBody(:,2) <= max(zFoil(:));
    
    if any(cond)
        
        ind(~cond) = [];

        % Upper and lower boundaries of radial body panels to include in plane
        % intersection calculation ie. potential panels that bridge point could
        % lie on. Wing not allowed to interfere with first or last radial
        % panels hence max(...,2), min(...,bDim-2) (Should be bDim-1?)
        % Pull in additional points just outside min (-1) and max (+1) z
        ub = max(min(ind)-1,2);
        lb = min(max(ind)+1,col-2);
        
    else
        % Upper/Lower boundaries may be entirely within two radial body
        % panels, thus no body points between min and max for cond. So
        % use all radial panels (bar first and last)
        ub = 2;
        lb = col-2;
        
    end
    
    j = ub:lb;
    
    interiorPoints = squeeze([foil(:,1,:);foil(:,end,:)]);
    exteriorPoints = squeeze([foil(:,2,:);foil(:,end-1,:)]);
    
    [fRow,fCol] = size(interiorPoints);
    bridge = zeros(fRow,fCol);
    
    innerBodyPoints = squeeze(bodyPoints(end,j,:));
    
    % Calculate bridge points ie. point on body where aerofoil joins
    for ii=1:fRow
        interiorPoint = interiorPoints(ii,:);
        exteriorPoint = exteriorPoints(ii,:);
        
        [inter,reason] = planeintersection(norm(j,:),innerBodyPoints,interiorPoint,exteriorPoint,xBody,yzBody,j);
        if isempty(inter)
            
            % Check what failed config looks like
%             aerofoil.Points = foil;
%             plotter({aftBody.Points,aerofoil.Points});
            
            success = false;
            return
        end
        bridge(ii,:) = inter;
    end
    
    if bool
        bridgeLower = reshape(bridge,[],1,3);
        bridgeUpper = zeros(fRow,1,3);
    else
        fRow = fRow/2;
        bridgeUpper = reshape(bridge(1:fRow,:),[],1,3);
        bridgeLower = reshape(bridge(fRow+1:end,:),[],1,3);
    end
    
    % Finds where other aerofoils interfere with said panels and futher discretises
    % only these panels
    zLE = bridge(1,3); % Leading edge z point
    
    zAbsDiff = abs(yzBody(:,2) - zLE);
    [~,upCut] = min(zAbsDiff);
    loCut = col - upCut;
    
    % Repmat not obviously necessary here, but needed to replace varying
    % number of elements in the body radial direction
    bridgeUpperMat = repmat(bridgeUpper,1,upCut);
    bridgeLowerMat = repmat(bridgeLower,1,loCut+1);
    
    rowUpper = 1:upCut;
    rowLower = (0:loCut) + upCut;
    
%     bodyUpper(:,:,1) = xBody = bodyPoints(:,1,1);
    
    %% Splitting body into separate upper and lower parts
    % Reshape back to 3D
    yzBodyUpper = reshape(yzBody(rowUpper,:),[],1,2);
    % Permute dimensions 1 & 2
    yzBodyUpper = permute(yzBodyUpper,[2,1,3]);
    % Repeat for number of upper bridge points
    yzBodyUpper = repmat(yzBodyUpper,fRow,1);
    % Add x coords to front from bridge
    bodyUpper(:,:,1) = bridgeUpperMat(:,:,1);
    bodyUpper(:,:,[2 3]) = yzBodyUpper;
    
    % Reshape back to 3D
    yzBodyLower = reshape(yzBody(rowLower,:),[],1,2);
    % Permute dimensions 1 & 2
    yzBodyLower = permute(yzBodyLower,[2,1,3]);
    % Repeat for number of lower bridge points
    yzBodyLower = repmat(yzBodyLower,fRow,1);
    % Add x coords to front from bridge
    bodyLower(:,:,1) = bridgeLowerMat(:,:,1);
    bodyLower(:,:,[2 3]) = yzBodyLower;
    
    % Find which body points to be replaced by upper/lower bridge
    conUpper = bodyUpper(:,:,3) < bridgeUpperMat(:,:,3);
    conLower = bodyLower(:,:,3) > bridgeLowerMat(:,:,3);
    
    % Set closest to wing bridges bodypoint lines equal to wing bridges
    conUpper(:,end) = true;
    conLower(:,1) = true;
    
    conUpper = repmat(conUpper,1,1,3);
    conLower = repmat(conLower,1,1,3);
    
    bodyUpper(conUpper) = bridgeUpperMat(conUpper);
    bodyLower(conLower) = bridgeLowerMat(conLower);
    
    %% Adding first and last chordwise body points to complete body
    
    bodyUpper = bodyUpper([1 1:end end],:,:);
    bodyLower = bodyLower([1 1:end end],:,:);
    
    % Only changes to be made are chordwise dimension (x)
    bodyUpper(1,:,1) = xBody(1);
    bodyLower(1,:,1) = xBody(1);
    
    bodyUpper(end,:,1) = xBody(end);
    bodyLower(end,:,1) = xBody(end);
    
    aftBody.Points = {bodyUpper,bodyLower};
    
    % Wing first and last point sets now replaced with upper and lower
    % bridge points found above
    foil(:,[1 end],:) = [bridgeUpper,bridgeLower];

    wingtail.Points = foil;
    
    %% CHECK: Does what it's supposed to?
    wetChord1 = diff(bridgeUpper([1 end],:,1));
    wetChord2 = wingtail.Chord(2);
    
    y1 = reshape(wingtail.Points(:,[1 end],2),[],1);
    y2 = reshape(wingtail.Points(:,[2 end-1],2),[],1);
    wetSpan = diff(mean([y1 y2],1));
    
    wetArea = ((wetChord2 + wetChord1)/2)*wetSpan;
    Taper = wetChord2/wetChord1;
    if Taper == inf
        Taper = 0;
    end
    wetMAC = (2/3)*wetChord1*((1+Taper+(Taper^2))/(1+Taper));
    
    wingtail.WetChord = [wetChord1, wingtail.Chord(2:end)];
    wingtail.WetArea = [wetArea, wingtail.Area(2:end)];
    wingtail.WetSpan = [wetSpan, wingtail.Span(2:end)];
    wingtail.WetMAC = [wetMAC, wingtail.MAC(2:end)];
    
    aerofoil(i) = wingtail;
end

success = true;