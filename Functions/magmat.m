function magnitudes = magmat(mat)
%% Calculates elemental magnitudes of 2D or 3D matrices

dims = size(mat);

if length(dims) == 3

    x = mat(:,:,1);
    y = mat(:,:,2);
    z = mat(:,:,3);
    
    magnitudes = (x.^2 + y.^2 + z.^2).^0.5;
    
elseif dims(2) == 3
    
    x = mat(:,1);
    y = mat(:,2);
    z = mat(:,3);
    
    magnitudes = (x.^2 + y.^2 + z.^2).^0.5;
    
else
    error('No setup for matrix of this form')
end
    