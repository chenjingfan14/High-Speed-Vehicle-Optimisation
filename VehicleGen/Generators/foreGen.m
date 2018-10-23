function fore = foreGen(foreProperties,Nose,aftbody,delta)

fore.Conical = true;
fore.Length = foreProperties(1);
fore.delta = delta;

xNose = Nose.x;
yNose = Nose.y;
zNose = Nose.z;

% Number of x-direction panels
xPanels = ceil(foreProperties(1)/delta);

fore.Points.Name = "forebody";

% Interpolate between nose points and start of aftbody
intx = [xNose; aftbody.x];
inty = [yNose; aftbody.y];
intz = [zNose; aftbody.z];

[dim,~] = size(xNose);

x_l = (0:1/xPanels:1)';
x = intx(dim,:) + (x_l).*(intx(end,:) - intx(dim,:));

if size(xNose,1) == 1
    % If nose is just a point, linear interpolate
    fore.Points.y = inty(1,:) + (x - x(1,:)).*((inty(2,:) - inty(1,:))./(intx(2,:)) - intx(1,:));
    fore.Points.z = intz(1,:) + (x - x(1,:)).*((intz(2,:) - intz(1,:))./(intx(2,:)) - intx(1,:));
    
else
    % If nose has non-zero parameters, use polynomial interpolation
    for j=1:size(inty,2)
        fore.Points.y(:,j) = interp1(intx(:,j),inty(:,j),x(:,j),'pchip');
        fore.Points.z(:,j) = interp1(intx(:,j),intz(:,j),x(:,j),'pchip');
    end
end

fore.Points.x = x;
fore.Points = xyztopoints(fore.Points);