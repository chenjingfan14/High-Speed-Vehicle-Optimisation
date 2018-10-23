function radial = radialpoints(struct)
%% Create single aftbody radial point matrix from all aftbody sections
% Since aftbody currently constant projection in x-direction, only first
% x-direction yz points are used since interpolation needs unique points

dim = length(struct);

xSecPoints = zeros(dim,1);
ySecPoints = zeros(dim,1);

% Find number of xPanels for each section of aftbody
% xSecPoints only needed if aftbody becomes variable in x-direction
% ySecPoints used to initialise radial point matrices
for i=1:dim
    aftSection = struct(i);
    [xSecPoints(i),ySecPoints(i)] = size(aftSection.x);
end

% Single point matrix size will equal sum of section sizes minus number of
% sections
yPoints = sum(ySecPoints) - (dim - 1);

[xAftbody,yAftbody,zAftbody] = deal(zeros(1,yPoints));
prevCol = 0;

for i=1:dim
    aftSection = struct(i);
    x = aftSection.x(1,:);
    y = aftSection.y(1,:);
    z = aftSection.z(1,:);
    
    colArray = 1:ySecPoints(i);
    
    % If a previous part has been analysed, use current parts first point 
    % and previous parts last point to average a single joining point
    if i > 1
        x(1) = (x(1) + xPrev)/2;
        y(1) = (y(1) + yPrev)/2;
        z(1) = (z(1) + zPrev)/2;
    end
    
    % If loop will continue, save final yz coords to average with first yz
    % coords of next part, and remove from current yz coords    
    if i ~= dim
        xPrev = x(end);
        yPrev = y(end);
        zPrev = z(end);
        
        x = x(1:end-1);
        y = y(1:end-1);
        z = z(1:end-1);
        
        colArray = colArray(1:end-1);
    end
    
    col = colArray + prevCol;
    prevCol = col(end);
    
    xAftbody(col) = x;
    yAftbody(col) = y;
    zAftbody(col) = z;
    
end

radial.x = xAftbody;
radial.y = yAftbody;
radial.z = zAftbody;