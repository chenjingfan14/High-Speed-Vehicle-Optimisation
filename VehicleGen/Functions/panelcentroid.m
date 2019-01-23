function centroidxyz = panelcentroid(x,y,z)
%% Calculate panel area
% Trapezoidal area found by spliting into two triangles and summing, works
% for any type of triangle

[row,col,dim3] = size(x);

row = 1:row - 1;

if nargin == 1
    
    col = 1:3:col - 3;
    
    if dim3 > 1
        % If 3D, ensure no singleton dimensions
        ax = squeeze(x(row, col, :));
        bx = squeeze(x(row + 1, col, :));
        cx = squeeze(x(row + 1, col + 3, :));
        dx = squeeze(x(row, col + 3, :));

        col = col + 1;

        ay = squeeze(x(row, col, :));
        by = squeeze(x(row + 1, col, :));
        cy = squeeze(x(row + 1, col + 3, :));
        dy = squeeze(x(row, col + 3, :));

        col = col + 1;

        az = squeeze(x(row, col, :));
        bz = squeeze(x(row + 1, col, :));
        cz = squeeze(x(row + 1, col + 3, :));
        dz = squeeze(x(row, col + 3, :));
    else
        
        ax = x(row, col);
        bx = x(row + 1, col);
        cx = x(row + 1, col + 3);
        dx = x(row, col + 3);

        col = col + 1;

        ay = x(row, col);
        by = x(row + 1, col);
        cy = x(row + 1, col + 3);
        dy = x(row, col + 3);

        col = col + 1;

        az = x(row, col);
        bz = x(row + 1, col);
        cz = x(row + 1, col + 3);
        dz = x(row, col + 3);
    end
    
elseif nargin == 3
    
    col = 1:col - 1;
    
    if dim3
        % If 3D, ensure no singleton dimensions
        ax = squeeze(x(row, col, :));
        bx = squeeze(x(row + 1, col, :));
        cx = squeeze(x(row + 1, col + 1, :));
        dx = squeeze(x(row, col + 1, :));

        ay = squeeze(y(row, col, :));
        by = squeeze(y(row + 1, col, :));
        cy = squeeze(y(row + 1, col + 1, :));
        dy = squeeze(y(row, col + 1, :));

        az = squeeze(z(row, col, :));
        bz = squeeze(z(row + 1, col, :));
        cz = squeeze(z(row + 1, col + 1, :));
        dz = squeeze(z(row, col + 1, :));
    else
        
        ax = x(row, col);
        bx = x(row + 1, col);
        cx = x(row + 1, col + 1);
        dx = x(row, col + 1);

        ay = x(row, col);
        by = x(row + 1, col);
        cy = x(row + 1, col + 1);
        dy = x(row, col + 1);

        az = x(row, col);
        bz = x(row + 1, col);
        cz = x(row + 1, col + 1);
        dz = x(row, col + 1);
    end

else
    
    error('Only one xyz matrix or three x, y, z matrices can be input')
end

centroid.x = (ax + bx + cx + dx)/4;
centroid.y = (ay + by + cy + dy)/4;
centroid.z = (az + bz + cz + dz)/4;

centroid = xyztopoints(centroid);

centroidxyz = centroid.xyz;