function radial = radialpoints(bodyCell)
%% Create single aftbody radial point matrix from all aftbody sections
% Since aftbody currently constant projection in x-direction, only first
% x-direction yz points are used since interpolation needs unique points

dim = length(bodyCell);

xSecPoints = zeros(dim,1);
ySecPoints = zeros(dim,1);

% Find number of xPanels for each section of aftbody
% xSecPoints only needed if aftbody becomes variable in x-direction
% ySecPoints used to initialise radial point matrices
for i=1:dim
    aftSection = bodyCell{i};
    [xSecPoints(i),ySecPoints(i),~] = size(aftSection);
end

% Single point matrix size will equal sum of section sizes minus number of
% sections
yPoints = sum(ySecPoints) - (dim - 1);

radial = zeros(yPoints,3);
prevRow = 0;

for i=1:dim
    
    % Squeezing to make things easier to index, expand back later
    xyz = squeeze(bodyCell{i}(1,:,:));
    rowArray = 1:ySecPoints(i);
    
    % If a previous part has been analysed, use current parts first point 
    % and previous parts last point to average a single joining point
    if i > 1
        
        xyz(1,:) = (xyz(1,:) + prev)/2;
    end
    
    % If loop will continue, save final yz coords to average with first yz
    % coords of next part, and remove from current yz coords    
    if i ~= dim
        
        prev = xyz(end,:);
        
        xyz = xyz(1:end-1,:);

        rowArray = rowArray(1:end-1);
    end
    
    col = rowArray + prevRow;
    prevRow = col(end);
    
    radial(col,:) = xyz;
    
end