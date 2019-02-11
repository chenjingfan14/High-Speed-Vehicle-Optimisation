function centroid = panelcentroid(xyz)
%% Calculate panel area
% Trapezoidal area found by spliting into two triangles and summing, works
% for any type of triangle

[row,col,~] = size(xyz);

row = 1:row - 1;
col = 1:col - 1;
    
a = xyz(row, col, :);
b = xyz(row + 1, col, :);
c = xyz(row + 1, col + 1, :);
d = xyz(row, col + 1, :);

centroid = (a + b + c + d)/4;