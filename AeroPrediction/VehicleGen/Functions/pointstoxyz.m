function new = pointstoxyz(old)
%% Translate single xyz point matrix to separate x, y, and z-coord matrices

new = old;
[row,col] = size(old);

for ii = 1:row*col
    
    xyz = old(ii).xyz;
    [~,dim] = size(xyz);
    
    X = 1:3:dim;
    Y = X + 1;
    Z = Y + 1;
    
    new(ii).x = xyz(:,X);
    new(ii).y = xyz(:,Y);
    new(ii).z = xyz(:,Z);
end