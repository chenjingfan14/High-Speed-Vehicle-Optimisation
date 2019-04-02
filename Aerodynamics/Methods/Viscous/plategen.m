function points = plategen()

%% Creating uniform/randomly perturbed panels
% Replace with size(points.x)
x = ones(2); % x,y points (x-1,y-1 panels)

% Creating x,y point matrices from 0 to [x y]-1
% Number of points
[dim1,dim2] = size(x);

% Number of panels
dim3 = dim1-1;
dim4 = dim2-1;

x = repmat((0:dim1-1)',1,dim2);
y = repmat((0:dim2-1),dim1,1);
z = zeros(dim1,dim2);

%% Random perturbation to make uniform grid chaotic, only inner points
randx = rand(dim3-1,dim4-1)*0.5 - 0.25;
randy = rand(dim3-1,dim4-1)*0.5 - 0.25;

x(2:end-1,2:end-1) = x(2:end-1,2:end-1) + randx;
y(2:end-1,2:end-1) = y(2:end-1,2:end-1) + randy;

%% Use for single quad
% randx = rand(2)*0.5 - 0.25;
% randy = rand(2)*0.5 - 0.25;
% 
% x = x + randx;
% y = y + randy;

%%

points(:,:,1) = x;
points(:,:,2) = y;
points(:,:,3) = z;
points = normals(points);
plotter(points)