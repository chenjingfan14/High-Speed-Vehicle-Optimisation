function vectors = crossmat(mat1,mat2)
%% Calculates elemental cross product x,y,z 3D vectors of x,y,z 3D matrices

[~,~,dim] = size(mat1);

if dim > 1
    
    x1 = mat1(:,:,1);
    y1 = mat1(:,:,2);
    z1 = mat1(:,:,3);
    
    x2 = mat2(:,:,1);
    y2 = mat2(:,:,2);
    z2 = mat2(:,:,3);
    
    vectors(:,:,1) = y1 .* z2 - z1 .* y2;
    vectors(:,:,2) = z1 .* x2 - x1 .* z2;
    vectors(:,:,3) = x1 .* y2 - y1 .* x2;
    
elseif numel(mat1) == 3
    
    x1 = mat1(1);
    y1 = mat1(2);
    z1 = mat1(3);
    
    x2 = mat2(1);
    y2 = mat2(2);
    z2 = mat2(3);
    
    vectors(1) = y1 .* z2 - z1 .* y2;
    vectors(2) = z1 .* x2 - x1 .* z2;
    vectors(3) = x1 .* y2 - y1 .* x2;
else
    error('Setup requried for this type of matrice input')
end