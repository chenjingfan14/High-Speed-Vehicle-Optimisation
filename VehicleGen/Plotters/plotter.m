function plotter(struct,varargin)

nCons = length(varargin);
conArray = 1:nCons;

plotDefs = repmat("",nCons,1);

for i = 1:nCons
    
    if isnumeric(varargin{i})
        
        plotDefs(i) = num2str(varargin{i});
    else
        plotDefs(i) = varargin(i);
    end
end

pos = figureposition();

figure

currentFig = gcf;

hold on
box on
set(currentFig, 'Position', pos)

[xlim, ylim, zlim] = deal([inf -inf]);

colour=[0.6 0.6 0.6];

if isstruct(struct) || iscell(struct)
    
    [~,nPoints] = size(struct);
else
    nPoints = 1;
end

for ii=1:nPoints
    
    if any(plotDefs == "triangle") && isfield(struct(ii),'Triangle')
            
        part = struct(ii).Triangle;
    else
        part = struct(ii);
    end
    
    if iscell(part)
        
        points = part{:};
        
    elseif isstruct(part)
        
        points = part.Points;
    else
        points = struct;
    end
    
    x = points(:,:,1);
    y = points(:,:,2);
    
    if size(points,3) == 2
        
        z = zeros(size(y));
    else
        z = points(:,:,3);
    end
    
    xlim = [min(xlim(1),min(x(:))) max(xlim(2),max(x(:)))];
    ylim = [min(ylim(1),min(y(:))) max(ylim(2),max(y(:)))];
    zlim = [min(zlim(1),min(z(:))) max(zlim(2),max(z(:)))];
    
    %%
    if any(plotDefs == "impact")
        
        x = part.Points(:,:,1);
        y = part.Points(:,:,2);
        z = part.Points(:,:,3);
        
        if isfield(part,'TriID')
            
            ID = part.TriID;
            
            x = x(ID);
            y = y(ID);
            z = z(ID);
            
        elseif isfield(part,'QuadID')
            
            ID = part.Quad;
            
            x = x(ID);
            y = y(ID);
            z = z(ID);
        end
        
        [row,~] = size(x);
        
        for i = 1:row
            
            if part.flow(i)
                
                colour = [1 0 0];
            else
                colour = [0.5 0.5 0.5];
            end
            fill3(x(i,:),y(i,:),z(i,:),colour);%,'EdgeColor','none')
        end
    elseif any(plotDefs == "wire")
        
        h = mesh(x,y,z,'edgecolor','k');
        
        % Removing legend entry for panels
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');    
        
    else
        
        if any(plotDefs == "triangle")
            
            [n,m] = size(x);
            
            if n > 1 && m == 3
                
                tri = zeros(n,m);
                dim = 1:numel(tri);
                
                tri(:) = dim;
                
            else
            
                mn = n * m;

                a = (1:mn - n)';

                t1 = [a, a+1, n+a];
                t2 = [a+1, n+a, n+a+1];

                delete = n:n:mn - n;

                t1(delete,:) = [];
                t2(delete,:) = [];

                tri = zeros(size(t1,1)*2,3);

                tri(1:2:end,:) = t1;
                tri(2:2:end,:) = t2;
            end
            
            tr = triangulation(tri,x(:),y(:),z(:));
            h = trisurf(tr,'LineWidth',0.05);
        else
            h = surf(x,y,z,'LineWidth',0.05);
        end
        if any(plotDefs == "nolines")
            
            set(h,'FaceColor',colour,'FaceLighting','flat','EdgeColor','none');
        else
            set(h,'FaceColor',colour,'FaceLighting','flat');
        end
        
        % Removing legend entry for panels
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end        
    
    %%
    if any(plotDefs == "centre")
        
        cx = part.centre(:,:,1);
        cy = part.centre(:,:,2);
        cz = part.centre(:,:,3);
        plot3(cx,cy,cz,'k.')
        
    end
    
    %%
    if any(plotDefs == "normal")
        
        cx = part.centre(:,:,1);
        cy = part.centre(:,:,2);
        cz = part.centre(:,:,3);
        
        nx = part.norm(:,:,1);
        ny = part.norm(:,:,2);
        nz = part.norm(:,:,3);
        plot3(nx + cx, ny + cy, nz + cz,'r.')
        
    end
    
    if any(plotDefs == "unit normal")
        
        cx = part.centre(:,:,1);
        cy = part.centre(:,:,2);
        cz = part.centre(:,:,3);
        
        % Dividing by constant factor to provide point just above surface
        nx = part.unitNorm(:,:,1)/10;
        ny = part.unitNorm(:,:,2)/10;
        nz = part.unitNorm(:,:,3)/10;
        plot3(nx + cx, ny + cy, nz + cz,'r.')
        
    end
    
    %%
    if any(plotDefs == "local")
        
        cx = part.centre(:,:,1);
        cy = part.centre(:,:,2);
        cz = part.centre(:,:,3);
        
        nx = part.unitNorm(:,:,1);
        ny = part.unitNorm(:,:,2);
        nz = part.unitNorm(:,:,3);
        
        unitTx = part.unitTang(:,:,1);
        unitTy = part.unitTang(:,:,2);
        unitTz = part.unitTang(:,:,3);
        
        unitSx = part.unitSurf(:,:,1);
        unitSy = part.unitSurf(:,:,2);
        unitSz = part.unitSurf(:,:,3);
        
        [row,col] = size(nx);
        
        for i = 1:row * col
            
            plot3([0 nx(i)/10]+cx(i),[0 ny(i)/10]+cy(i),[0 nz(i)/10]+cz(i),'r')
            plot3([0 unitTx(i)/10]+cx(i),[0 unitTy(i)/10]+cy(i),[0 unitTz(i)/10]+cz(i),'b')
            plot3([0 unitSx(i)/10]+cx(i),[0 unitSy(i)/10]+cy(i),[0 unitSz(i)/10]+cz(i),'g')
        end
        
    end
    
    %%
    if any(plotDefs == "CoP")
        
        if isfield(points,'CoP')
            CoP = points.CoP;
        else
            
            which = plotDefs == "CoP";
            which = conArray(which);
            
            CoP = str2double(plotDefs(which+1));
        end
        
        plot3(CoP(1),CoP(2),CoP(3),'r.','markers',72)
    end
    
    %%
    if any(plotDefs == "double")
        if any(plotDefs == "wire")
            
            mesh(x,-y,z,'edgecolor','k');
            
        else
            
            h = surf(x,-y,z,'LineWidth',0.05);
            
            if any(plotDefs == "nolines")
                set(h,'FaceColor',colour,'FaceLighting','flat','EdgeColor','none');
            else
                set(h,'FaceColor',colour,'FaceLighting','flat');
            end
        end
        
        ylim = [min(ylim(1),min(-y(:))) max(ylim(2),max(y(:)))];
    end
    
    if any(plotDefs == "title")
        
        which = plotDefs == "title";
        which = conArray(which);
        
        title(plotDefs{which+1});
    end
    
end

if xlim(1) == xlim(2)
    
    xlim(2) = xlim(2) + 1;
end

if ylim(1) == ylim(2)
    
    ylim(2) = ylim(2) + 1;
end

if zlim(1) == zlim(2)
    
    zlim(2) = zlim(2) + 1;
end

axis('equal',[xlim ylim zlim])
xlabel('x, m','Interpreter','latex','FontSize',14);
ylabel('y, m','Interpreter','latex','FontSize',14);
zlabel('z, m','Interpreter','latex','FontSize',14);
set(gca,'TickLabelInterpreter','latex','FontSize', 14);
hold off

if any(plotDefs == "pause")
    disp('Simulation paused until key is pressed')
    pause
end