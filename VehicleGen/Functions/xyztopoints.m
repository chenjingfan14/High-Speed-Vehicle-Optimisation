function points = xyztopoints(x,y,z)

if size(x,2) == 1
    
    x = repmat(x,1,size(y,2));
end

if size(y,1) == 1
    
    y = repmat(y,size(x,1),1);
end
if size(z,1) == 1
    
    z = repmat(z,size(x,1),1);
end

points(:,:,1) = x;
points(:,:,2) = y;
points(:,:,3) = z;