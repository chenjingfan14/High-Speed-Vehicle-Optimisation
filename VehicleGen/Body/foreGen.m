function fore = foreGen(foreProperties,nosePoints,aftbody,delta)

fore.Name = "forebody";
fore.Conical = true;
fore.Length = foreProperties(1);
fore.delta = delta;

% Number of x-direction panels
xPanels = ceil(foreProperties(1)/delta);

% fore.Points.Name = "forebody";

% Interpolate between nose points and start of aftbody
int = [nosePoints; reshape(aftbody,1,[],3)];

[row,col,~] = size(nosePoints);

x_l = (0:1/xPanels:1)';
x = int(row,:,1) + (x_l).*(int(end,:,1) - int(row,:,1));

if row == 1
    % If nose is just a point, linear interpolate
    yz = int(1,:,[2 3]) + (x - x(1,:)).*((int(2,:,[2 3]) - int(1,:,[2 3]))./(x(end,:)) - x(1,:));
    
else
    % If nose has non-zero parameters, use polynomial interpolation
    for j = col:-1:1
        
        yz(:,j,:) = interp1(int(:,j,1),int(:,j,[2 3]),x(:,j),'pchip');
    end
end

fore.Points(:,:,1) = x;
fore.Points(:,:,[2 3]) = yz;