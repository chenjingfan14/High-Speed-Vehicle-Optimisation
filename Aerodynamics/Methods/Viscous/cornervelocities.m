function partCell = cornervelocities(partCell,flow)

% clear all
% close all
%
% flow = flowparameters();
% points = plategen();

Vinf = flow.U;

for ii = length(partCell):-1:1
    
    time = 0;
    
    part = partCell{ii};
    
    norm = part.norm;
    tang = part.unitTang;
    surf = part.unitSurf;
    L = part.L;
    
    %% PUT OUTSIDE OF FUNCTION
    % Non-dimensionalise so that all parts lie in positive x,y & z
    xyz = part.Points - min(min(part.Points));
    part = normals(xyz);
    
    xyz = part.Points;
    
    centre = part.centre;
    panelArea = part.area;
    
    %%
    % Creating x,y point matrices from 0 to [x y]-1
    % Number of points
    [dim1,dim2,~] = size(xyz);
    
    % Number of panels
    dim3 = dim1 - 1;
    dim4 = dim2 - 1;
    
    %% Centre point velocities
    
    T = crossmat(norm, permute(Vinf,[3 1 2]));
    Vcentre = crossmat(T, norm);
    
    % Translate to corner velocity
    Vcorner = L * Vcentre(:);
    
    Vc = reshape(Vcorner,size(part.Points));
    
    %% Initialise streamline tracing
    projxyz = xyz(:,:,1);
    projxyz(:,:,2) = (xyz(:,:,2).^2 + xyz(:,:,3).^2).^0.5;
    
    projCentre = centre(:,:,1);
    projCentre(:,:,2) = (centre(:,:,2).^2 + centre(:,:,3).^2).^0.5;
    
    ex = [1 0 0];
    ey = [0 1 0];
    ez = [0 0 1];
    
    ew = reshape(norm,[],3);
    eu = reshape(tang,[],3);
    ev = reshape(surf,[],3);
    
    %% THIS IS FOR POINTS NOT CENTRES IE +1
        
    inject = projCentre(1,:,:);
    area = panelArea(1,:);
    
    delete = area == 0;
    inject(:,delete,:) = [];
    
    cx = centre(:,:,1);
    cy = centre(:,:,2);
    cz = centre(:,:,3);
    
    % What is this. Fix.
    ou = cx - cx;
    ov = cy - projCentre(:,:,2);
    ow = cz - cz;
    
    plotter(projxyz)
    
    %%
    
    Ax = xyz(1:end-1,1:end-1,1);
    Bx = xyz(2:end,1:end-1,1);
    Cx = xyz(2:end,2:end,1);
    Dx = xyz(1:end-1,2:end,1);
    
    Ay = projxyz(1:end-1,1:end-1,2);
    By = projxyz(2:end,1:end-1,2);
    Cy = projxyz(2:end,2:end,2);
    Dy = projxyz(1:end-1,2:end,2);
    
    %%
    rows = ones(1,dim4);
    cols = 1:dim4;
    
    rows(delete) = [];
    cols(delete) = [];
    
    nRows = dim3*100;
    
    % Initial matrix ID, first leading edge panels
    ID = rows + (cols - 1) * dim3;
    [~,tdim,~] = size(inject);
    tcols = 1:tdim;
    timeID = 2 + (tcols-1)*nRows;
    
    [xp,yp] = deal(nan(nRows,tdim));
    
    x1 = Ax(ID);
    x2 = Bx(ID);
    x3 = Cx(ID);
    x4 = Dx(ID);
    
    y1 = Ay(ID);
    y2 = By(ID);
    y3 = Cy(ID);
    y4 = Dy(ID);
    
    c1 = [x1;y2];
    
    %% Velocity definition
    
    % Translate to panel corner velocities
    Vx = Vc(:,:,1);
    Vy = Vc(:,:,2);
    Vz = Vc(:,:,3);
    
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
    
    timestep = ones(1,tdim)*0.0001;
    
    V1 = [Vx1(ID); Vy1(ID)];
    V2 = [Vx2(ID); Vy2(ID)];
    V3 = [Vx3(ID); Vy3(ID)];
    V4 = [Vx4(ID); Vy4(ID)];
    
    inject = permute(inject,[3 2 1]);
    
    xp(1,:) = inject(1,:);
    yp(1,:) = inject(2,:);
    
    %% FOR LOOP THIS SHIT
    
    for count = 1:tdim
    
        IDi = ID(count);
    
        eui = eu(IDi,:);
        evi = ev(IDi,:);
        ewi = ew(IDi,:);
    
        R = eui.*ex' + evi.*ey' + ewi.*ez';
    
        Rv = R*[0 0 0]';
    
        pNew(1,count) = Rv(1) + ou(IDi) + cx(IDi);
        ypnew(1,count) = Rv(2) + ov(IDi) + cy(IDi);
        zpnew(1,count) = Rv(3) + ow(IDi) + cz(IDi);
    
        pNew(1,count) = Rv(1) + cx(IDi);
        ypnew(1,count) = Rv(2) + cy(IDi);
        zpnew(1,count) = Rv(3) + cz(IDi);
    
    end
    %%
    
    figure(1)
    hold on
    axis equal
    plot(xp(1,:),yp(1,:),'r*')
    hold off
    
    plotter(part,"centre")
    
    figure(2)
    hold on
    axis equal
    plot3(pNew,ypnew,zpnew,'r*')
    hold off
    
    Vp1 = V1 + (inject - c1).*((V4 - V1)./(y4 - y1));
    
    con1 = Vp1 == inf;
    con2 = Vp1 == -inf;
    con3 = isnan(Vp1);
    
    Vp1(con1 | con2 | con3) = 1;
    
    mini = -0.00001;
    maxi = 1.00001;
    
    xp1 = xp(1,:);
    yp1 = yp(1,:);
    xp2 = xp1 + 1e-16;
    yp2 = yp1 + 1e-16;
    
    [xCross,yCross,xPrevCross,yPrevCross,test] = deal(nan(4,tdim));
    
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
    
    nx = norm(:,:,1);
    ny = norm(:,:,2);
    nz = norm(:,:,3);
    
    nx0 = nx(ID);
    ny0 = ny(ID);
    nz0 = nz(ID);
    
    realtime = zeros(1,tdim);
    model = 'barycentricmodel';
    particleArray = 1:tdim;
    
    while any(~stopCon)
        
        pt = [xp1;yp1];
        
        % RK4 integration
        k1 = timestep.*feval(model, [pt; Vp1] ,timestep,x1,x2,x3,x4,y1,y2,y3,y4,V1,V2,V3,V4);
        k2 = timestep.*feval(model, [pt; Vp1] + k1/2 ,timestep,x1,x2,x3,x4,y1,y2,y3,y4,V1,V2,V3,V4);
        k3 = timestep.*feval(model, [pt; Vp1] + k2/2 ,timestep,x1,x2,x3,x4,y1,y2,y3,y4,V1,V2,V3,V4);
        k4 = timestep.*feval(model, [pt; Vp1] + k3 ,timestep,x1,x2,x3,x4,y1,y2,y3,y4,V1,V2,V3,V4);
        
        delta = (k1 + 2*k2 + 2*k3 + k4)/6;
        
        deltax = delta(1,:);
        deltay = delta(2,:);
        
        Vp2 = Vp1 + delta(3:4,:);
        
        xpt = xp1 + deltax;
        ypt = yp1 + deltay;
        
        xp2 = xpt;
        yp2 = ypt;
        
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
        test(:,ydec) = [x1(ydec);y1(ydec);x2(ydec);y2(ydec)];
        
        % Particle cannot cross back over line it has just crossed, must be
        % different otherwise no intersection
        ydec = ydec & any(test ~= yPrevCross,1);
        
        yCross(:,ydec) = test(:,ydec);
        
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
        
        test(:,xinc) = [x2(xinc);y2(xinc);x3(xinc);y3(xinc)];
        xinc = xinc & any(test ~= xPrevCross,1);
        xCross(:,xinc) = test(:,xinc);
        
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
        
        test(:,yinc) = [x4(yinc);y4(yinc);x3(yinc);y3(yinc)];
        yinc = yinc & any(test ~= yPrevCross,1);
        yCross(:,yinc) = test(:,yinc);
        
        %% Has particle crossed between panel points 4 & 1
        Px4 = ((xp1.*yp2 - yp1.*xp2).*(x4-x1) - (xp1-xp2).*(x4.*y1 - y4.*x1))./...
            ((xp1-xp2).*(y4-y1) - (yp1-yp2).*(x4-x1));
        
        Py4 = ((xp1.*yp2 - yp1.*xp2).*(y4-y1) - (yp1-yp2).*(x4.*y1 - y4.*x1))./...
            ((xp1-xp2).*(y4-y1) - (yp1-yp2).*(x4-x1));
        
        t4 = ((xp1-x4).*(y4-y1) - (yp1-y4).*(x4-x1))./...
            ((xp1-xp2).*(y4-y1) - (yp1-yp2).*(x4-x1));
        
        u4 = -((xp1-xp2).*(yp1-y4) - (yp1-yp2).*(xp1-x4))./...
            ((xp1-xp2).*(y4-y1) - (yp1-yp2).*(x4-x1));
        
        xdec = t4 >= mini & t4 <= maxi & u4 >= mini & u4 <= maxi;
        
        test(:,xdec) = [x1(xdec);y1(xdec);x4(xdec);y4(xdec)];
        xdec = xdec & any(test ~= xPrevCross,1);
        
        test(:,xdec) = [x4(xdec);y4(xdec);x1(xdec);y1(xdec)];
        xdec = xdec & any(test ~= xPrevCross,1);
        xCross(:,xdec) = test(:,xdec);
        
        if any(xdec)
           
            xdec
        end
        
        %% Updates for next timestep
        % Line crossed in current timestep becomes previous crossing point for
        % future timesteps
        xPrevCross(:,xdec) = xCross(:,xdec);
        xPrevCross(:,xinc) = xCross(:,xinc);
        yPrevCross(:,ydec) = yCross(:,ydec);
        yPrevCross(:,yinc) = yCross(:,yinc);
        
        % Updates panel that particle is within
        xdec1 = xdec & rows ~= 1;
        xdec2 = xdec & rows == 1;
        
        if any(xdec1)
            
            rows(xdec1) = rows(xdec1) - 1;
        end
        
        if any(xdec2)
            
            cols(xdec2) = dim4 - cols(xdec2) + 1;
        end
        
        if any(ydec)
            
            ydec
        end
        
        if any(xinc)
            
            xinc
        end
        
        if any(yinc)
            
            yinc
        end
        
        rows(xinc) = rows(xinc) + 1;
        cols(yinc) = cols(yinc) + 1;
        cols(ydec) = cols(ydec) - 1;
        
        ID = rows + (cols-1)*dim3;
        
        % Starting point for next timestep is intersection point. Removes error
        % of particle crossing current panel AND a t+1 panel line which will
        % not be picked up by simulation. However time and velocity must now be
        % updated for such situations
        xpt(ydec) = Px1(ydec);
        xpt(xinc) = Px2(xinc);
        xpt(yinc) = Px3(yinc);
        xpt(xdec) = Px4(xdec);

        ypt(ydec) = Py1(ydec);
        ypt(xinc) = Py2(xinc);
        ypt(yinc) = Py3(yinc);
        ypt(xdec) = Py4(xdec);
        
        % Streamline points updated
        xp(timeID) = xpt;
        yp(timeID) = ypt;
        
%         for count = 1:sum(~stopCon)
%             
%             IDi = ID(count);
%             timeIDi = timeID(count);
%             
%             eui = eu(IDi,:);
%             evi = ev(IDi,:);
%             ewi = ew(IDi,:);
%             
%             R = eui.*ex' + evi.*ey' + ewi.*ez';
%             
% %             Rv = R*[xpt(count)-cx(IDi) ypt(count)-yzCentre(IDi) 0]';
%             Rv = R*[0 0 0]';
%             
%             pNew(timeIDi) = Rv(1) + ou(IDi);
%             ypnew(timeIDi) = Rv(2) + ov(IDi);
%             zpnew(timeIDi) = Rv(3) + ow(IDi);
%             
%             pNew(timeIDi) = Rv(1) + cx(IDi) + ou(IDi);
%             ypnew(timeIDi) = Rv(2) + cy(IDi) + ov(IDi);
%             zpnew(timeIDi) = Rv(3) + cz(IDi) + ow(IDi);
%             
%             pNew(timeIDi) = Rv(1) + cx(IDi);
%             ypnew(timeIDi) = Rv(2) + cy(IDi);
%             zpnew(timeIDi) = Rv(3) + cz(IDi);
%             
%         end
        
        %% CHECK/DELETE
        %         euID = [eux(ID) euy(ID) euz(ID)];
        %         evID = [evx(ID) evy(ID) evz(ID)];
        %         ewID = [ewx(ID) ewy(ID) ewz(ID)];
        %
        %         R = euID*ex' + evID*ey' + ewID*ez';
        %
        %         Rv = R*[xp; yp; zeros(1,sum(~stopCon))];
        %
        %         xpnew(timeID) = Rv(1,:) + ou;
        %         ypnew(timeID) = Rv(2,:) + ov;
        %         zpnew(timeID) = Rv(3,:) + ow;
        
        %% Plotting every timestep
        figure(1)
        hold on
        for j=1:sum(~stopCon)
            dim = timeID(j)-1:timeID(j);
            plot(xp(dim),yp(dim),'r')
            %         plot(xp(dim),yp(dim),'r*')
        end
        hold off
        
%         figure(2)
%         hold on
%         for j=1:sum(~stopCon)
%             dim = timeID(j)-1:timeID(j);
%             plot3(pNew(dim),ypnew(dim),zpnew(dim),'r*')
%         end
%         hold off
        
        %%
        cross = xinc | yinc | ydec;
        
        pdiff = ((xpt-xp1).^2 + (ypt-yp1).^2).^0.5;
        
        Vpint = Vp1 + (xpt-xp1).*((Vp2-Vp1)./(xp2-xp1));
        
        Vp2(:,cross) = Vpint(:,cross);
        
        xp1 = xpt;
        yp1 = ypt;
        Vp1 = Vp2;
        
        Vavg = (Vp1+Vp2)./2;
        Vmag = (Vavg(1,:).^2 + Vavg(2,:).^2).^0.5;
        
        realtime(particleArray) = realtime(particleArray) + pdiff./Vmag;
        
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
            
            test = test(:,~stopCon);
            
            xp1 = xp1(:,~stopCon);
            yp1 = yp1(:,~stopCon);
            Vp1 = Vp1(:,~stopCon);
            
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
        zeroArea = panelArea(ID) == 0;
        
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
            
            test = test(:,~stopCon);
            
            xp1 = xp1(:,~stopCon);
            yp1 = yp1(:,~stopCon);
            Vp1 = Vp1(:,~stopCon);
            
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
        
        % Unit vectors so magnitudes = 1
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
            
            test = test(:,~stopCon);
            
            xp1 = xp1(:,~stopCon);
            yp1 = yp1(:,~stopCon);
            Vp1 = Vp1(:,~stopCon);
            
            cross = cross(~stopCon);
            
            timestep = timestep(~stopCon);
            timeID = timeID(~stopCon);
        end
        %%
        
        nx0 = nx1;
        ny0 = ny1;
        nz0 = nz1;
        
        % New corner points and velocities for next step
        x1 = Ax(ID);
        x2 = Bx(ID);
        x3 = Cx(ID);
        x4 = Dx(ID);
        
        y1 = Ay(ID);
        y2 = By(ID);
        y3 = Cy(ID);
        y4 = Dy(ID);
        
        V1 = [Vx1(ID); Vy1(ID)];
        V4 = [Vx4(ID); Vy4(ID)];
        V2 = [Vx2(ID); Vy2(ID)];
        V3 = [Vx3(ID); Vy3(ID)];
        
%         swap23 = y2 - y3 > 0;
%         xtemp = x2(swap23);
%         ytemp = y2(swap23);
%         Vtemp = V2(:,swap23);
%         
%         x2(swap23) = x3(swap23);
%         x3(swap23) = xtemp;
%         
%         y2(swap23) = y3(swap23);
%         y3(swap23) = ytemp;
%         
%         V2(:,swap23) = V3(:,swap23);
%         V3(:,swap23) = Vtemp;
%         
%         swap14 = y1 - y4 > 0;
%         xtemp = x1(swap14);
%         ytemp = y1(swap14);
%         Vtemp = V1(:,swap14);
%         
%         x1(swap14) = x4(swap14);
%         x4(swap14) = xtemp;
%         
%         y1(swap14) = y4(swap14);
%         y4(swap14) = ytemp;
%         
%         V1(:,swap14) = V4(:,swap14);
%         V4(:,swap14) = Vtemp;
        
        timeID = timeID + 1;
        
        %% Variable timestep update
        %         if any(cross)
        %
        %             xc = [x1;x2;x3;x4];
        %             yc = [y1;y2;y3;y4];
        %
        %             xPrev = xPrevCross(:,cross);
        %             yPrev = yPrevCross(:,cross);
        %
        %             timestep(cross) = estimateexit(xp1(cross),yp1(cross),Vp1(:,cross),xc(:,cross),yc(:,cross),xPrev,yPrev,mini,maxi);
        %
        %         end
        %
        %         timestep(~cross) = timestep(~cross)*2;
        
        time = time + max(timestep);
        
        if time >= 0.01
            
            break
        end
        
    end
    
    close all
    
%     figure(1)
%     hold on
%     plot(xp,yp,'r')
%     hold off
    
%     figure(2)
%     hold on
%     plot3(pNew,ypnew,zpnew,'r')
%     hold off
    
%     partCell{ii} = part;
    
end