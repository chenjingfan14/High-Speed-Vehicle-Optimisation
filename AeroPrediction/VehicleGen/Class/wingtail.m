classdef wingtail
    
    properties
        Conical = false
        delta = 0.1
        Offset
        xPanels = 50
        % "Linear", "Cosine" or "HalfCosine"
        Distribution = "Cosine"
        Points
        Chord
        Area
        Span
        MAC
        WetChord
        WetArea
        WetSpan
        WetMAC
        Partitions
        Boolean
        
    end
    
    methods
        function obj=wingtail(dihedral,semispan,chord,sweep,offset,sections)
            
            % If any aerofoil section span = 0, delete all corresponding partition
            % properties
            dim = length(semispan);
            cond = semispan == 0;
            semispan(cond) = [];
            sweep(cond) = [];
            
            ind = (1:dim)+1;
            ind = ind(cond);
            chord(ind) = [];
            
            for j = ind
                sections{j} = '';
            end
            
            sections(strcmp('',sections)) = [];
            
            %%
            
            xPanels = obj.xPanels;
            X = (0:xPanels)';
            c = chord;
            partitions = sum(semispan>0);
            
            [bh,Area,cbar] = deal(zeros(1,partitions));
            xLE = zeros(1,partitions+1);
            
            % Is aerofoil vertical tail
            di = dihedral*pi/180;
            boolean = di == pi/2;
            
            Lam = sweep*pi/180;
            picker = obj.Distribution;
            
            prev=0;
            for i=1:partitions
                bh(i) = prev+semispan(i);
                xLE(i+1) = xLE(i)+semispan(i)*tan(Lam(i));
                Area(i) = ((c(i+1)+c(i))/2)*semispan(i);
                Taper = c(i+1)/c(i);
                if Taper == inf
                    Taper = 0;
                end
                cbar(i) = (2/3)*c(i)*((1+Taper+(Taper^2))/(1+Taper));
                
                prev = bh(i);
                
            end
            
            obj.Area = Area;
            obj.WetArea = Area;
            obj.MAC = cbar;
            obj.WetMAC = cbar;
            
            y = [0,bh]*cos(di);
            z = offset(2)+([0,bh]*sin(di));
            
            %% Wing x-distribution
            switch picker
                case "Linear"
                    xnorm = X/max(X);
                case "Cosine"
                    xnorm = 0.5*(1-cos((X*pi)/max(X)));
                case "HalfCosine"
                    xnorm = 1-cos((X*(pi/2))/max(X));
            end
            
            %% Aerofoil Section x-Discretisation
            
            [FoilUp,FoilLo] = deal(zeros(size(xnorm,1),partitions+1));

            for i=1:length(sections)
                
                Sec = sections{i};
                
                dim = size(Sec,1);
                half = round((dim+1)/2);
                                
                upper = unique(Sec(1:half,:),'stable','rows');
                lower = unique(Sec(half:dim,:),'stable','rows');
                
                % Interpolate z given section xz coords and defined 
                % x-discretisation 
                FoilUp(:,i) = interp1(upper(:,1),upper(:,2),xnorm,'pchip'); % interpolates sample data wrt defined chorwise points
                FoilLo(:,i) = interp1(lower(:,1),lower(:,2),xnorm,'pchip');
            end
            
            % Calculates all xpoints with leading edge sweep offset
            [x_u,x_l] = deal(xnorm*c + xLE);
            
            if di < pi/4 && di > -pi/4
                y_u = y+FoilUp*-sin(di).*c;
                z_u = z+FoilUp*cos(di).*c;
                y_l = y+FoilLo*-sin(di).*c;
                z_l = z+FoilLo*cos(di).*c;
                
                % Interp/Extrapolate to align root with centreline
                if di~=0
                    z_u(:,1) = z_u(:,1)+(-y_u(:,1).*((z_u(:,2)-z_u(:,1))./(y_u(:,2)-y_u(:,1))));
                    y_u(:,1) = 0;
                    z_l(:,1) = z_l(:,1)+(-y_l(:,1).*((z_l(:,2)-z_l(:,1))./(y_l(:,2)-y_l(:,1))));
                    y_l(:,1) = 0;
                end
            else
                % This doesn't work
            end
            
            if ~isequal(z_u(end,:),z_l(end,:))
                x_u(end+1,:) = x_l(end,:);
                y_u(end+1,:) = y_l(end,:);
                z_u(end+1,:) = z_l(end,:);
                x_l(end+1,:) = x_l(end,:);
                y_l(end+1,:) = y_l(end,:);
                z_l(end+1,:) = z_l(end,:);
            end
            
            if boolean
                wing.x = x_l;
                wing.y = y_l;
                wing.z = z_l;
            else
                wing.x = [x_u, fliplr(x_l)];
                wing.y = [y_u, fliplr(y_l)];
                wing.z = [z_u, fliplr(z_l)];
            end
            
            obj.Chord = chord;
            obj.WetChord = chord;
            obj.Span = semispan;
            obj.WetSpan = semispan;
            obj.Offset = offset;
            obj.Points = xyztopoints(wing);
            obj.Points.Name = "aerofoil";
            obj.Partitions = partitions;
            obj.Boolean = boolean;
            
            %%
%             plotter(obj.Points);
        end
    end
    
end

