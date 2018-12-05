function particleviewer(assemblyProperties,num)

prev = 0;
dimArray = zeros(100,1);
assemblyPoints = struct('x', [], 'y', [], 'z', [], 'xyz', []);
for i=1:length(assemblyProperties)
    dim = (1:length(assemblyProperties{i}.Points))+prev;
    assemblyPoints(dim) = [assemblyProperties{i}.Points];
    assemblyProperties{i}.Points = [];
    dimArray(i) = dim(end)-prev;
    prev = dim(end);
end

% new = rotate(assemblyPoints,MAC,flow.alpha);
% rotatedPoints = normals(new,flow);

title = ['Configuration ' num2str(num)];

plotter(assemblyPoints,"title",title);