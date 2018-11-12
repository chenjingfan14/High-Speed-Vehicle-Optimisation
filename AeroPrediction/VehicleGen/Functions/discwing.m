function new = discwing(old)

new = old;
target = old.delta;
dim1 = length(new);

for ii=1:dim1
    
    wingtail = old(ii);
    
    bool = wingtail.Boolean;
    x = wingtail.Points.x;
    y = wingtail.Points.y;
    z = wingtail.Points.z;
    
    [rows,c1] = size(x);
    
    ybar = mean(y,1);
    zbar = mean(z,1);
    
    yzbar = (ybar.^2 + zbar.^2).^0.5;
    
    if bool
        % fix for tail parts
    else
        
        % Base spanwise discretisation on either upper or lower aerofoil,
        % pick based on smallest of span of the two
        upperDiff = abs(yzbar(c1/2) - yzbar(1));
        lowerDiff = abs(yzbar((c1/2)+1) - yzbar(end));
        
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
        
    end
    
    % Points = Total number of panels + 1
    cols = sum(numPanel)+1;
    
    [xp,yp,zp] = deal(zeros(rows,cols));
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
        yi = y0 + (y1 - y0).*(panelArray/panels);
        
        z0 = z(:,j);
        z1 = z(:,j+1);
        zi = z0 + (z1 - z0).*(panelArray/panels);
        
        %% x-coords based on combination of y and z
        
        yz0 = (y0.^2 + z0.^2).^0.5;
        yz1 = (y1.^2 + z1.^2).^0.5;
        yzi = yz0 + (yz1 - yz0).*(panelArray/panels);
        
        xi = x(:,j) + (yzi - yz0).*((x(:,j+1) - x(:,j))./(yz1 - yz0));
        
        cols = panelArray + dim;
        
        xp(:,cols) = xi;
        yp(:,cols) = yi;
        zp(:,cols) = zi;
        
        dim = cols(end);
        count = count + 1;
    end
    
    %% Update points to include discretisation
    
    new(ii).Points.x = xp;
    new(ii).Points.y = yp;
    new(ii).Points.z = zp;
    new(ii).Points = xyztopoints(new(ii).Points);
end

%% Plot undiscretised wing and new discretised version
% plotter(old.Points)
% plotter(new.Points)

end