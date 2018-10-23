function intstreamline(points,~)

[~,dim] = size(points);

for i=1:dim
    
    part = points(i);
    
    %% Creating uniform/randomly perturbed panels
    
%     xyz = part.xyz;
%     offset = xyz(1,1:3);
    
    x = part.x; % x,y points (x-1,y-1 panels)
    y = part.y;% - offset(2);
    z = part.z;% - offset(3);
    
    theta = atan2(y,z);
    
        at = theta(1:end-1,1:end-1);
    bt = theta(2:end,1:end-1);
    ct = theta(2:end,2:end);
    dt = theta(1:end-1,2:end);
    
    centrex = part.centre(:,1:3:end);
    centrey = part.centre(:,2:3:end);
    centrez = part.centre(:,3:3:end);
    
    centreyz = (centrey.^2 + centrez.^2).^0.5;
    theta = atan2(centrey,centrez);
    
    area = part.area;
    
    % Creating x,y point matrices from 0 to [x y]-1
    % Number of points
    [dim1,dim2] = size(x);
    
    % Number of panels
    dim3 = dim1-1;
    dim4 = dim2-1;
    
    norm = part.unitNorm;
    
    X = 1:3:dim4*3;
    Y = X + 1;
    Z = Y + 1;
    
    nx = norm(:,X);
    ny = norm(:,Y);
    nz = norm(:,Z);
    
    nTransMag = (ny.^2 + nz.^2).^0.5;
    
    nyTrans = ny./nTransMag;
    nzTrans = nz./nTransMag;
    
    theta = atan2(nyTrans,nzTrans);
    
    test = (nyTrans.^2 + nzTrans.^2).^0.5;
    
    xinject = centrex(1,:);
    yinject = centrey(1,:);
    zinject = centrez(1,:);
    theta1 = atan2(yinject,zinject);
    
    %% Transform yz coords
    yz = (y.^2 + z.^2).^0.5;
    
    transpoints = part;
    transpoints.y = yz;
    transpoints.z = zeros(dim1,dim2);
    transpoints = xyztopoints(transpoints);
    
    plotter(transpoints)
    
    %%
    
    ax = x(1:end-1,1:end-1);
    bx = x(2:end,1:end-1);
    cx = x(2:end,2:end);
    dx = x(1:end-1,2:end);
    
    ay = yz(1:end-1,1:end-1);
    by = yz(2:end,1:end-1);
    cy = yz(2:end,2:end);
    dy = yz(1:end-1,2:end);
    
    % Edge mid points
    %     chordx = (x(:,2:end) + x(:,1:end-1))/2;
    %     chordy = (yz(:,2:end) + yz(:,1:end-1))/2;
    %
    %     spanx = (x(2:end,:) + x(1:end-1,:))/2;
    %     spany = (yz(2:end,:) + yz(1:end-1,:))/2;
    
    % Initial injection points
    %     xinject = (chordx(1,:) + chordx(2,:))/2;
    %     yinject = (chordy(1,:) + chordy(2,:))/2;
    
    yinject = (yinject.^2 + zinject.^2).^0.5;
    
    %%
    rows = ones(1,dim4);
    cols = 1:dim4;
    
    rows = [rows]; %1:dim3];
    cols = [cols]; %ones(1,dim3)];
    
    nRows = dim3*100;
    
    % Initial matrix ID, first leading edge panels
    ID = rows + (cols-1)*dim3;
    [~,tdim] = size(xinject);
    tcols = 1:tdim;
    timeID = 2 + (tcols-1)*nRows;
    
    dummy = nan(nRows,tdim);
    
    xp = dummy;
    yp = dummy;
    
    xpnew = dummy;
    ypnew = dummy;
    zpnew = dummy;
    
    thetaTrans = dummy;
    thetaTrans(1,:) = theta1;
    
    x1 = ax(ID);
    x2 = bx(ID);
    x3 = cx(ID);
    x4 = dx(ID);
    
    y1 = ay(ID);
    y2 = by(ID);
    y3 = cy(ID);
    y4 = dy(ID);
    
    the1 = at(ID);
    the2 = bt(ID);
    the3 = ct(ID);
    the4 = dt(ID);
    
    c1 = [x1;y2];
    
    %% Velocity definition
    
    Vc = part.Vc;
    
    % Translate to panel corner velocities
    Vx = Vc(:,1:3:end);
    Vy = Vc(:,2:3:end);
    Vz = Vc(:,3:3:end);
    
%     Vx(:) = 0;
%     Vy(:) = 0;
%     Vz(:) = 0;
    
    %% Transform yz velocity
    Vyneg = Vy < 0;
    Vzneg = Vz < 0;
    
    Vysquared = Vy.^2;
    Vzsquared = Vz.^2;
    
    Vysquared(Vyneg) = -Vysquared(Vyneg);
    Vzsquared(Vzneg) = -Vzsquared(Vzneg);
    
    Vy = (Vysquared + Vzsquared).^0.5;
    
    Re = real(Vy) > 0;
    Im = imag(Vy) > 0;
    
    % New variable?
    Vy(Re) = real(Vy(Re));
    Vy(Im) = -imag(Vy(Im));
    
    Vx1 = Vx(1:end-1,1:end-1);
    Vx2 = Vx(2:end,1:end-1);
    Vx3 = Vx(2:end,2:end);
    Vx4 = Vx(1:end-1,2:end);
    
    Vy1 = Vy(1:end-1,1:end-1);
    Vy2 = Vy(2:end,1:end-1);
    Vy3 = Vy(2:end,2:end);
    Vy4 = Vy(1:end-1,2:end);
    
    %% Initialise integration
    
    stopCon = false(1,tdim);
    
    timestep = ones(1,tdim)*0.00001;
    
    V1 = [Vx1(ID); Vy1(ID)];
    V2 = [Vx2(ID); Vy2(ID)];
    V3 = [Vx3(ID); Vy3(ID)];
    V4 = [Vx4(ID); Vy4(ID)];
    
    xp(1,:) = xinject;
    yp(1,:) = yinject;
    
    xpnew(1,:) = xinject;
    ypnew(1,:) = yinject.*sin(theta1);% + offset(2);
    zpnew(1,:) = yinject.*cos(theta1);% + offset(3);
    
    figure(1)
    hold on
    axis equal
    plot(xp(1,:),yp(1,:),'r*')
    hold off
    
    plotter(points)
    
    figure(2)
    hold on
    axis equal
    plot3(xpnew(1,:),ypnew(1,:),zpnew(1,:),'r*')
    hold off
    
    %     Vp1(1,:) = ones(1,dim4);
    %     Vp1(2,:) = -dim4/2:(dim4-1)/2;
    Vp1 = V1+(yinject-c1).*((V4-V1)./(y4-y1));
    Vp1 = zeros(2,tdim);
    
    con1 = Vp1 == inf;
    con2 = Vp1 == -inf;
    con3 = isnan(Vp1);
    
    Vp1(con1 | con2 | con3) = 1;
    
    mini = -0.00001;
    maxi = 1.00001;
    
    xp1 = xpnew(1,:);
    yp1 = ypnew(1,:);
    xp2 = xpnew(1,:) + 1e-16;
    yp2 = ypnew(1,:) + 1e-16;
    
    xCross = nan(4,tdim);
    yCross = nan(4,tdim);
    xPrevCross = nan(4,tdim);
    yPrevCross = nan(4,tdim);
    
    %% 2D
    %% Has particle crossed between panel points 1 & 2
    % Defines whether intersection point is between two given line segments
    t1 = ((xp1-x1).*(y1-y2) - (yp1-y1).*(x1-x2))./...
        ((xp1-xp2).*(y1-y2) - (yp1-yp2).*(x1-x2));
    
    u1 = -((xp1-xp2).*(yp1-y1) - (yp1-yp2).*(xp1-x1))./...
        ((xp1-xp2).*(y1-y2) - (yp1-yp2).*(x1-x2));
    
    % Is intersection is between both line segments
    ydec = t1 >= mini & t1 <= maxi & u1 >= mini & u1 <= maxi;
    
    % Which panel line is being crossed
    yCross(:,ydec) = [x1(ydec);y1(ydec);x2(ydec);y2(ydec)];
    
    % Particle cannot cross back over line it has just crossed, must be
    % different otherwise no intersection
    ydec = ydec & any(yCross ~= yPrevCross,1);
    
    %% Has particle crossed between panel points 2 & 3
    t2 = ((xp1-x2).*(y2-y3) - (yp1-y2).*(x2-x3))./...
        ((xp1-xp2).*(y2-y3) - (yp1-yp2).*(x2-x3));
    
    u2 = -((xp1-xp2).*(yp1-y2) - (yp1-yp2).*(xp1-x2))./...
        ((xp1-xp2).*(y2-y3) - (yp1-yp2).*(x2-x3));
    
    xinc = t2 >= mini & t2 <= maxi & u2 >= mini & u2 <= maxi;
    xCross(:,xinc) = [x2(xinc);y2(xinc);x3(xinc);y3(xinc)];
    xinc = xinc & any(xCross ~= xPrevCross,1);
    
    %% Has particle crossed between panel points 3 & 4
    t3 = ((xp1-x3).*(y3-y4) - (yp1-y3).*(x3-x4))./...
        ((xp1-xp2).*(y3-y4) - (yp1-yp2).*(x3-x4));
    
    u3 = -((xp1-xp2).*(yp1-y3) - (yp1-yp2).*(xp1-x3))./...
        ((xp1-xp2).*(y3-y4) - (yp1-yp2).*(x3-x4));
    
    yinc = t3 >= mini & t3 <= maxi & u3 >= mini & u3 <= maxi;
    yCross(:,yinc) = [x4(yinc);y4(yinc);x3(yinc);y3(yinc)];
    yinc = yinc & any(yCross ~= yPrevCross,1);
    
    %% Updates for next timestep
    % Line crossed in current timestep becomes previous crossing point for
    % future timesteps
    xPrevCross(:,xinc) = xCross(:,xinc);
    yPrevCross(:,ydec) = yCross(:,ydec);
    yPrevCross(:,yinc) = yCross(:,yinc);
    
    nx0 = nx(ID);
    ny0 = ny(ID);
    nz0 = nz(ID);
    
    realtime = zeros(1,tdim);
    model = 'barycentricmodel';
    particleArray = 1:tdim;
    integrate = true;
    
    if integrate
        while any(~stopCon)
            
            pt = [xp1;yp1];
            
            % RK4 integration
            k1 = timestep.*feval(model, [pt; Vp1; theta1] ,timestep,x1,x2,x3,x4,y1,y2,y3,y4,V1,V2,V3,V4,the1,the2,the3,the4);
            k2 = timestep.*feval(model, [pt; Vp1; theta1] + k1/2 ,timestep,x1,x2,x3,x4,y1,y2,y3,y4,V1,V2,V3,V4,the1,the2,the3,the4);
            k3 = timestep.*feval(model, [pt; Vp1; theta1] + k2/2 ,timestep,x1,x2,x3,x4,y1,y2,y3,y4,V1,V2,V3,V4,the1,the2,the3,the4);
            k4 = timestep.*feval(model, [pt; Vp1; theta1] + k3 ,timestep,x1,x2,x3,x4,y1,y2,y3,y4,V1,V2,V3,V4,the1,the2,the3,the4);
            
            delta = (k1 + 2*k2 + 2*k3 + k4)/6;
            
            xp2 = xp1 + delta(1,:);
            yp2 = yp1 + delta(2,:);
            Vp2 = Vp1 + delta(3:4,:);
            theta2 = theta1 + delta(5,:);
            
            %% 2D
            %% Has particle crossed between panel points 1 & 2
            % Intersection point of two infinite lines
            Px1 = ((xp1.*yp2 - yp1.*xp2).*(x1-x2) - (xp1-xp2).*(x1.*y2 - y1.*x2))./...
                ((xp1-xp2).*(y1-y2) - (yp1-yp2).*(x1-x2));
            
            Py1 = ((xp1.*yp2 - yp1.*xp2).*(y1-y2) - (yp1-yp2).*(x1.*y2 - y1.*x2))./...
                ((xp1-xp2).*(y1-y2) - (yp1-yp2).*(x1-x2));
            
            % Defines whether intersection point is between two given line segments
            t1 = ((xp1-x1).*(y1-y2) - (yp1-y1).*(x1-x2))./...
                ((xp1-xp2).*(y1-y2) - (yp1-yp2).*(x1-x2));
            
            u1 = -((xp1-xp2).*(yp1-y1) - (yp1-yp2).*(xp1-x1))./...
                ((xp1-xp2).*(y1-y2) - (yp1-yp2).*(x1-x2));
            
            % Is intersection is between both line segments
            ydec = t1 >= mini & t1 <= maxi & u1 >= mini & u1 <= maxi;
            
            % Which panel line is being crossed
            yCross(:,ydec) = [x1(ydec);y1(ydec);x2(ydec);y2(ydec)];
            
            % Particle cannot cross back over line it has just crossed, must be
            % different otherwise no intersection
            ydec = ydec & any(yCross ~= yPrevCross,1);
            
            %% Has particle crossed between panel points 2 & 3
            Px2 = ((xp1.*yp2 - yp1.*xp2).*(x2-x3) - (xp1-xp2).*(x2.*y3 - y2.*x3))./...
                ((xp1-xp2).*(y2-y3) - (yp1-yp2).*(x2-x3));
            
            Py2 = ((xp1.*yp2 - yp1.*xp2).*(y2-y3) - (yp1-yp2).*(x2.*y3 - y2.*x3))./...
                ((xp1-xp2).*(y2-y3) - (yp1-yp2).*(x2-x3));
            
            t2 = ((xp1-x2).*(y2-y3) - (yp1-y2).*(x2-x3))./...
                ((xp1-xp2).*(y2-y3) - (yp1-yp2).*(x2-x3));
            
            u2 = -((xp1-xp2).*(yp1-y2) - (yp1-yp2).*(xp1-x2))./...
                ((xp1-xp2).*(y2-y3) - (yp1-yp2).*(x2-x3));
            
            xinc = t2 >= mini & t2 <= maxi & u2 >= mini & u2 <= maxi;
            xCross(:,xinc) = [x2(xinc);y2(xinc);x3(xinc);y3(xinc)];
            xinc = xinc & any(xCross ~= xPrevCross,1);
            
            %% Has particle crossed between panel points 3 & 4
            Px3 = ((xp1.*yp2 - yp1.*xp2).*(x3-x4) - (xp1-xp2).*(x3.*y4 - y3.*x4))./...
                ((xp1-xp2).*(y3-y4) - (yp1-yp2).*(x3-x4));
            
            Py3 = ((xp1.*yp2 - yp1.*xp2).*(y3-y4) - (yp1-yp2).*(x3.*y4 - y3.*x4))./...
                ((xp1-xp2).*(y3-y4) - (yp1-yp2).*(x3-x4));
            
            t3 = ((xp1-x3).*(y3-y4) - (yp1-y3).*(x3-x4))./...
                ((xp1-xp2).*(y3-y4) - (yp1-yp2).*(x3-x4));
            
            u3 = -((xp1-xp2).*(yp1-y3) - (yp1-yp2).*(xp1-x3))./...
                ((xp1-xp2).*(y3-y4) - (yp1-yp2).*(x3-x4));
            
            yinc = t3 >= mini & t3 <= maxi & u3 >= mini & u3 <= maxi;
            yCross(:,yinc) = [x4(yinc);y4(yinc);x3(yinc);y3(yinc)];
            yinc = yinc & any(yCross ~= yPrevCross,1);
            
            %% Updates for next timestep
            % Line crossed in current timestep becomes previous crossing point for
            % future timesteps
            xPrevCross(:,xinc) = xCross(:,xinc);
            yPrevCross(:,ydec) = yCross(:,ydec);
            yPrevCross(:,yinc) = yCross(:,yinc);
            
            % Updates panel that particle is within
            rows(xinc) = rows(xinc) + 1;
            cols(yinc) = cols(yinc) + 1;
            cols(ydec) = cols(ydec) - 1;
            
            ID = rows + (cols-1)*dim3;
            
            % Starting point for next timestep is intersection point. Removes error
            % of particle crossing current panel AND a t+1 panel line which will
            % not be picked up by simulation. However time and velocity must now be
            % updated for such situations
            xpt = xp2;
            ypt = yp2;
            
            xpt(ydec) = Px1(ydec);
            xpt(xinc) = Px2(xinc);
            xpt(yinc) = Px3(yinc);
            
            ypt(ydec) = Py1(ydec);
            ypt(xinc) = Py2(xinc);
            ypt(yinc) = Py3(yinc);
            
            %%
            cross = xinc | yinc | ydec;
            
            pt1 = ((xpt-xp1).^2 + (ypt-yp1).^2).^0.5;
            p21 = ((xp2-xp1).^2 + (yp2-yp1).^2).^0.5;
            
%             Vpint = Vp1 + (xpt-xp1).*((Vp2-Vp1)./(xp2-xp1));
%             thetaint = theta1 + (xpt-xp1).*((theta2-theta1)./(xp2-xp1));
            
            Vpint = Vp1 + pt1.*((Vp2-Vp1)./p21);
            thetaint = theta1 + pt1.*((theta2-theta1)./p21);
            thetaint = theta1 + (ypt-yp1).*((theta2-theta1)./(yp2-yp1));
            
            Vp2(:,cross) = Vpint(:,cross);
            theta2(cross) = thetaint(:,cross);
            theta2 = thetaint;
            
            xp1 = xpt;
            yp1 = ypt;
            Vp1 = Vp2;
            theta1 = theta2;
            
            Vavg = (Vp1+Vp2)./2;
            Vmag = (Vavg(1,:).^2 + Vavg(2,:).^2).^0.5;
            
%             thetaID = barycentric([xpt;ypt],x1,x2,x3,x4,y1,y2,y3,y4,the1,the2,the3,the4);
            thetaTrans(timeID) = theta1;
            centrey(ID);
            centrez(ID);
            
            ypnew(timeID) = (ypt-yp1).*sin(theta(ID)) + centrey(ID);
            zpnew(timeID) = (ypt-yp1).*cos(theta(ID)) + centrez(ID);
            
            % Streamline points updated
            xp(timeID) = xpt;
            yp(timeID) = ypt;
            
            xpnew(timeID) = xpt;
            ypnew(timeID) = centrey(ID) + ypt.*sin(theta(ID));
            zpnew(timeID) = centrez(ID) + ypt.*cos(theta(ID));
            
            %% Plotting every timestep
            figure(1)
            hold on
            for j=1:sum(~stopCon)
                dim = timeID(j)-1:timeID(j);
                plot(xp(dim),yp(dim),'r')
%                 plot(xp(dim),yp(dim),'r*')
            end
            hold off
            
            figure(2)
            hold on
            for j=1:sum(~stopCon)
                dim = timeID(j)-1:timeID(j);
                plot3(xpnew(dim),ypnew(dim),zpnew(dim),'r')
%                 plot3(xpnew(dim),ypnew(dim),zpnew(dim),'r*')
            end
            hold off
            
            realtime(particleArray) = realtime(particleArray) + pt1./Vmag;
            
            %% Stop particle
            % Stopping conditions: When particle reaches surface boundary
            maxRow = rows >= dim1;
            maxCol = cols >= dim2;
            minCol = cols <= 0;
            
            % Stop particle if any of above conditions met
            stopCon = maxRow | maxCol | minCol;
            
            %% FIRST UPDATE
            % Remove rows & cols of particles that have reached surface
            % boundary
            if any(stopCon)
                rows = rows(~stopCon);
                cols = cols(~stopCon);
                
                xinc = xinc(~stopCon);
                yinc = yinc(~stopCon);
                ydec = ydec(~stopCon);
                
                ID = ID(~stopCon);
                
                particleArray = particleArray(~stopCon);
                
                xCross = xCross(:,~stopCon);
                yCross = yCross(:,~stopCon);
                
                xPrevCross = xPrevCross(:,~stopCon);
                yPrevCross = yPrevCross(:,~stopCon);
                
                xp1 = xp1(~stopCon);
                yp1 = yp1(~stopCon);
                Vp1 = Vp1(:,~stopCon);
                theta1 = theta1(~stopCon);
                
                nx0 = nx0(~stopCon);
                ny0 = ny0(~stopCon);
                nz0 = nz0(~stopCon);
                
                cross = cross(~stopCon);
                
                timestep = timestep(~stopCon);
                timeID = timeID(~stopCon);
            end
            %% Zero area next cell condition
            % If turning angle too large then assume line has detached and
            % stop
            
            zeroArea = area(ID) == 0;
            
            % Updates panel that particle is within
            rows(xinc & zeroArea) = rows(xinc & zeroArea) + 1;
            cols(yinc & zeroArea) = cols(yinc & zeroArea) + 1;
            cols(ydec & zeroArea) = cols(ydec & zeroArea) - 1;
            
            % Stopping conditions: When particle reaches surface boundary
            maxRow = rows >= dim1;
            maxCol = cols >= dim2;
            minCol = cols <= 0;
            
            % Stop particle if any of above conditions met
            stopCon = maxRow | maxCol | minCol;
            
            %% SECOND UPDATE
            % Remove rows & cols of particles that have reached surface
            % boundary
            if any(stopCon)
                rows = rows(~stopCon);
                cols = cols(~stopCon);
                
                ID = ID(~stopCon);
                
                particleArray = particleArray(~stopCon);
                
                xCross = xCross(:,~stopCon);
                yCross = yCross(:,~stopCon);
                
                xPrevCross = xPrevCross(:,~stopCon);
                yPrevCross = yPrevCross(:,~stopCon);
                
                xp1 = xp1(~stopCon);
                yp1 = yp1(~stopCon);
                Vp1 = Vp1(:,~stopCon);
                theta1 = theta1(~stopCon);
                
                nx0 = nx0(~stopCon);
                ny0 = ny0(~stopCon);
                nz0 = nz0(~stopCon);
                
                cross = cross(~stopCon);
                
                timestep = timestep(~stopCon);
                timeID = timeID(~stopCon);
            end
            
            %% Turning angle stopping condition
            % Calculate dot product and magnitudes to work out turning
            % angle
            
            nx1 = nx(ID);
            ny1 = ny(ID);
            nz1 = nz(ID);
            
            dot = nx0.*nx1 + ny0.*ny1 + nz0.*nz1;
            
            % Unit vectors so magnitude = 1
%             mag0 = (nx0.^2 + ny0.^2 + nz0.^2).^0.5;
%             mag1 = (nx1.^2 + ny1.^2 + nz1.^2).^0.5;
            
            turnAngle = real(acos(dot));
            
            % Fix with actual angle
            stopCon = turnAngle > 60.*pi/180;
            
            %% THIRD UPDATE
            if any(stopCon)
                % Remove rows & cols of particles that have reached surface
                % boundary
                rows = rows(~stopCon);
                cols = cols(~stopCon);
                
                ID = ID(~stopCon);
                
                nx1 = nx(ID);
                ny1 = ny(ID);
                nz1 = nz(ID);
                
                particleArray = particleArray(~stopCon);
                
                xCross = xCross(:,~stopCon);
                yCross = yCross(:,~stopCon);
                
                xPrevCross = xPrevCross(:,~stopCon);
                yPrevCross = yPrevCross(:,~stopCon);
                
                xp1 = xp1(~stopCon);
                yp1 = yp1(~stopCon);
                Vp1 = Vp1(:,~stopCon);
                theta1 = theta1(~stopCon);
                
                cross = cross(~stopCon);
                
                timestep = timestep(~stopCon);
                timeID = timeID(~stopCon);
            end
            %%
            
            nx0 = nx1;
            ny0 = ny1;
            nz0 = nz1;
            
            % New corner points and velocities for next step
            x1 = ax(ID);
            x2 = bx(ID);
            x3 = cx(ID);
            x4 = dx(ID);
            
            y1 = ay(ID);
            y2 = by(ID);
            y3 = cy(ID);
            y4 = dy(ID);
            
            the1 = at(ID);
            the2 = bt(ID);
            the3 = ct(ID);
            the4 = dt(ID);
            
            V1 = [Vx1(ID); Vy1(ID)];
            V2 = [Vx2(ID); Vy2(ID)];
            V3 = [Vx3(ID); Vy3(ID)];
            V4 = [Vx4(ID); Vy4(ID)];
            
            timeID = timeID + 1;
            
            %% Variable timestep update
            if any(cross)
                
                
%             
%                 xc = [x1;x2;x3;x4];
%                 yc = [y1;y2;y3;y4];
%                 
%                 xPrev = xPrevCross(:,cross);
%                 yPrev = yPrevCross(:,cross);
%                 
%                 timestep(cross) = estimateexit(xp1(cross),yp1(cross),Vp1(:,cross),xc(:,cross),yc(:,cross),xPrev,yPrev,mini,maxi);
% 
            end
%             
%             timestep(~cross) = timestep(~cross)*2;
            
        end
    else
        
    end
end
figure(2)
hold on
plot3(xpnew,ypnew,zpnew,'r')
% plot3(xpnew,ypnew,zpnew,'r*')
hold off