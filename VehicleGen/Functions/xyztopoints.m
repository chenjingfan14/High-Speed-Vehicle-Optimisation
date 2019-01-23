function new = xyztopoints(old)

[row,col] = size(old);
new = old;

for ii=1:row*col
    
    x = old(ii).x;
    y = old(ii).y;
    z = old(ii).z;
    
    if size(x,2) == 1
        x = repmat(x,1,size(y,2));
        new(ii).x = x;
    end
    if size(y,1) == 1
        y = repmat(y,size(x,1),1);
        new(ii).y = y;
    end 
    if size(z,1) == 1
        z = repmat(z,size(x,1),1);
        new(ii).z = z;
    end
    
    [dim1,dim2,dim3] = size(x);
    
    new(ii).xyz = zeros(dim1,dim2*3,dim3);
    
    X = 1:3:dim2*3;
    Y = X + 1;
    Z = Y + 1;
    
    new(ii).xyz(:,X,:) = x;
    new(ii).xyz(:,Y,:) = y;
    new(ii).xyz(:,Z,:) = z;

end