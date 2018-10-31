function plotter(pointStruct,varargin)

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

h = findobj('type','figure');
num = length(h)+1;

num6 = ceil(num/6);

if num6 > 1
    num = num - 6*(num6-1);
end

n = num./(1:6);
isWhole = rem(n,1) == 0;
lowestWhole = max(n(isWhole));

switch lowestWhole
    case 1
        pos = [0, 560, 640, 440];
    case 3 
        pos = [640, 560, 640, 440];
    case 5
        pos = [1280, 560, 640, 440];
    case 2
        pos = [0, 40, 640, 440];
    case 4
        pos = [640, 40, 640, 440];
    case 6
        pos = [1280, 40, 640, 440];
        
%     case 1
%         pos = [0, -150, 640, 440];
%     case 3 
%         pos = [640, -150, 640, 440];
%     case 5
%         pos = [1280, -150, 640, 440];
%     case 2
%         pos = [0, 380, 640, 440];
%     case 4
%         pos = [640, 380, 640, 440];
%     case 6
%         pos = [1280, 380, 640, 440];
        
end

figure
hold on
set(gcf, 'Position', pos)
axis equal
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')

% r=[1,0,0];
% g=[0,1,0];
% b=[0,0,1];

colour=[0.6 0.6 0.6];

[~,nPoints] = size(pointStruct);

for ii=1:nPoints
    
    points = pointStruct(ii);
    
    x = points.x;
    y = points.y;
    z = points.z;
    
    [~,dim2] = size(x);
    dim3 = (dim2-1)*3;
    
    X = 1:3:dim3;
    Y = X + 1;
    Z = Y + 1;
    
    %%
    if any(plotDefs == "impact")
        
        centre = points.centre;
        xyz = points.xyz;
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
        
        cx = points.centre(:,X);
        cy = points.centre(:,Y);
        cz = points.centre(:,Z);
        plot3(cx,cy,cz,'k*')
        
    end
    
    %%
    if any(plotDefs == "normals")
        
        cx = points.centre(:,X);
        cy = points.centre(:,Y);
        cz = points.centre(:,Z);
        
        nx = points.norm(:,X);
        ny = points.norm(:,Y);
        nz = points.norm(:,Z);
        plot3(nx + cx, ny + cy, nz + cz,'r*')
        
    end
    
    %%
    if any(plotDefs == "local")
        
        cx = points.centre(:,X);
        cy = points.centre(:,Y);
        cz = points.centre(:,Z);
        
        nx = points.norm(:,X);
        ny = points.norm(:,Y);
        nz = points.norm(:,Z);
        
        unitTx = points.unitTang(:,X);
        unitTy = points.unitTang(:,Y);
        unitTz = points.unitTang(:,Z);
        
        unitSx = points.unitSurf(:,X);
        unitSy = points.unitSurf(:,Y);
        unitSz = points.unitSurf(:,Z);
        
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
            
            CoP = str2num(plotDefs(which+1));
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