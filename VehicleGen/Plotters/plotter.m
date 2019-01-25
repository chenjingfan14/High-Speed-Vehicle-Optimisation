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
figNum = currentFig.Number;

hold on
set(currentFig, 'Position', pos)
axis equal
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')

% r=[1,0,0];
% g=[0,1,0];
% b=[0,0,1];

colour=[0.6 0.6 0.6];

[~,nPoints] = size(struct);

for ii=1:nPoints
    
    part = struct(ii);
    
    if iscell(part)
        
        points = part{:};
    else
        points = part.Points;
    end
    
    x = points(:,:,1);
    y = points(:,:,2);
    z = points(:,:,3);
    
    [~,dim2] = size(points);
    dim3 = dim2 - 3;
    
    X = 1:3:dim3;
    Y = X + 1;
    Z = Y + 1;
    
    %%
    if any(plotDefs == "impact")
        
        centre = part.centre;
        xyz = part.xyz;
        [x,y] = size(centre);
        
        for i=1:x
            for j=1:3:y
                
                p = [xyz(i:i+1,j:j+2); xyz(i+1:-1:i,j+3:j+5)];
                
                if points.flow(i,(j+2)/3) == 1
                    colour = [1 0 0];
                else
                    colour = [0.5 0.5 0.5];
                end
                fill3(p(:,1),p(:,2),p(:,3),colour);%,'EdgeColor','none');
            end
        end
    elseif any(plotDefs == "wire")
        
        mesh(x,y,z,'edgecolor','k');
        
    else
        
        h = surf(x,y,z,'LineWidth',0.05);
        
        if any(plotDefs == "nolines")
            set(h,'FaceColor',colour,'FaceLighting','flat','EdgeColor','none');
        else
            set(h,'FaceColor',colour,'FaceLighting','flat');
        end
    end
    
    %%
    if any(plotDefs == "centre")
        
        cx = part.centre(:,:,1);
        cy = part.centre(:,:,2);
        cz = part.centre(:,:,3);
        plot3(cx,cy,cz,'k*')
        
    end
    
    %%
    if any(plotDefs == "normals")
        
        cx = part.centre(:,:,1);
        cy = part.centre(:,:,2);
        cz = part.centre(:,:,3);
        
        nx = part.norm(:,:,1);
        ny = part.norm(:,:,2);
        nz = part.norm(:,:,3);
        plot3(nx + cx, ny + cy, nz + cz,'r*')
        
    end
    
    if any(plotDefs == "unit normals")
        
        cx = part.centre(:,:,1);
        cy = part.centre(:,:,2);
        cz = part.centre(:,:,3);
        
        % Dividing by constant factor to provide point just above surface
        nx = part.unitNorm(:,:,1)/10;
        ny = part.unitNorm(:,:,2)/10;
        nz = part.unitNorm(:,:,3)/10;
        plot3(nx + cx, ny + cy, nz + cz,'r*')
        
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
        
        for i = 1:row*col
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
    end
    
    if any(plotDefs == "title")
        
        which = plotDefs == "title";
        which = conArray(which);
        
        title(plotDefs{which+1});
    end
    
end
hold off

if any(plotDefs == "pause")
    disp('Simulation paused until key is pressed')
    pause
end