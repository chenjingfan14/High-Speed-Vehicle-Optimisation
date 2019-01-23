function new = pointstoxyz(old)
%% Translate single xyz point matrix to separate x, y, and z-coord matrices

new = old;
[row,col,~] = size(old);

for ii = 1:row*col
    
    xyz = old(ii).xyz;
    [~,dim,dim3] = size(xyz);
    
    X = 1:3:dim;
    Y = X + 1;
    Z = Y + 1;
    
    new(ii).x = squeeze(xyz(:,X,:));
    new(ii).y = xyz(:,Y,:);
    new(ii).z = xyz(:,Z,:);
    
    if dim3 > 1
        % Ensure no singleton dimensions
        new(ii).x = squeeze(xyz(:,X,:));
        new(ii).y = squeeze(xyz(:,Y,:));
        new(ii).z = squeeze(xyz(:,Z,:));
    else
        
        new(ii).x = xyz(:,X,:);
        new(ii).y = xyz(:,Y,:);
        new(ii).z = xyz(:,Z,:);
    end
end