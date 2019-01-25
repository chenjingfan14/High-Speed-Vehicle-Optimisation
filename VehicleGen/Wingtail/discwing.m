function wingtail = discwing(wingtail)

target = wingtail.delta;
dim1 = length(wingtail);
controlSurf = wingtail.Control.Surf;

for ii=1:dim1
    
    bool = wingtail(ii).Boolean;
    xyz = wingtail(ii).Points;
    
    x = xyz(:,:,1);
    y = xyz(:,:,2);
    z = xyz(:,:,3);
    
    [rows,c1] = size(x);
    
    ybar = mean(y,1);
    zbar = mean(z,1);
    
    yzbar = (ybar.^2 + zbar.^2).^0.5;
    
    if bool
        % fix for tail parts
    else
        
        % Base spanwise discretisation on either upper or lower aerofoil,
        % pick based on smallest of span of the two
        upperDiff = abs(diff(yzbar([1 c1/2])));
        lowerDiff = abs(diff(yzbar([end (c1/2)+1])));
        
        minDiff = min(upperDiff,lowerDiff);
        
        if upperDiff == minDiff % Base on upper
            use = yzbar(1:c1/2);
            real = ceil(abs(use(2:end) - use(1:end-1))/target);
            numPanel = [real, fliplr(real)];
        else % Base on lower
            use = yzbar((c1/2)+1:end);
            real = ceil(abs(use(2:end) - use(1:end-1))/target);
            numPanel = [fliplr(real), real];
        end
        
        numPanel(numPanel == 0) = 1;
        
    end
    
    % Points = Total number of panels + 1
    col = sum(numPanel)+1;
    
    newControl = false(1,(col+1)/2);
    
    [xp,yp,zp] = deal(zeros(rows,col));
    [dim,count] = deal(1);
    
    for j=1:c1-1
        
        % Skip wing tip, no need to discretise here
        if ~bool && j == c1/2
            dim = dim+1;
            continue
        end
        
        panels = numPanel(count);
        panelArray = 0:panels;
        
        %% Find intermediate yz-coords
        y0 = y(:,j);
        y1 = y(:,j+1);
        
        if all(y0 == y1)
            yi = repmat(y0,1,panels+1);
        else
            yi = y0 + (y1 - y0).*(panelArray/panels);
        end
        
        z0 = z(:,j);
        z1 = z(:,j+1);
        
        if all(z0 == z1)
            zi = repmat(z0,1,panels+1);
        else
            zi = z0 + (z1 - z0).*(panelArray/panels);
        end
        
        %% x-coords based on combination of y and z
        
        x0 = x(:,j);
        x1 = x(:,j+1);
        
        if all(x0 == x1)
            xi = repmat(x0,1,panels+1);
        else
            yz0 = (y0.^2 + z0.^2).^0.5;
            yz1 = (y1.^2 + z1.^2).^0.5;
            yzi = yz0 + (yz1 - yz0).*(panelArray/panels);
            xi = x(:,j) + (yzi - yz0).*((x(:,j+1) - x(:,j))./(yz1 - yz0));
        end
        
        if  j < c1/2 && controlSurf(count)
            newControl(cols) = true;
        end
        
        if panels == 0
            cols = dim+1;
        else
            cols = panelArray + dim;
        end
        
        
        xp(:,cols) = xi;
        yp(:,cols) = yi;
        zp(:,cols) = zi;
        
        dim = cols(end);
        count = count + 1;
    end
    
    first = find(newControl,1,'first');
    newControl(first) = false;

    %% Update points to include discretisation
    wingtail(ii).Control.Surf = [newControl, fliplr(newControl)];
    
    wingtail(ii).Points = [];
    wingtail(ii).Points(:,:,1) = xp;
    wingtail(ii).Points(:,:,2) = yp;
    wingtail(ii).Points(:,:,3) = zp;
end

%% Plot undiscretised wing and new discretised version
% plotter(old.Points)
% plotter(new.Points)

end