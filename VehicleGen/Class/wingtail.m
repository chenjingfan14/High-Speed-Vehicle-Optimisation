classdef wingtail
    
    properties
        Name
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
        Dihedral
        MAC
        WetChord
        WetArea
        WetSpan
        WetMAC
        Partitions
        Control
        Boolean
        
    end
    
    methods
        function obj=wingtail(dihedral,semispan,chord,sweep,sections,control)
            
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
            nParts = sum(semispan > 0);
            nSecs = nParts + 1;
            
            [Area,cbar] = deal(zeros(1,nParts));
            [bh,xLE] = deal(zeros(1,nParts+1));
            
            % Is aerofoil vertical tail
            di = dihedral*pi/180;
            boolean = di == pi/2;
            
            Lam = sweep*pi/180;
            picker = obj.Distribution;
            
            for i=1:nParts
                bh(i+1) = bh(i) + semispan(i);
                xLE(i+1) = xLE(i) + semispan(i)*tan(Lam(i));
                Area(i) = 0.5*(chord(i+1) + chord(i)) * semispan(i);
                Taper = chord(i+1)/chord(i);
                if Taper == inf
                    Taper = 0;
                end
                cbar(i) = (2/3)*chord(i)*((1 + Taper + (Taper^2))/(1 + Taper));
                
            end
            
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
            
            [FoilUp,FoilLo] = deal(zeros(size(xnorm,1),nParts+1));

            for i=1:nSecs
                
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
            
            %% Control surfaces
            if nargin == 7
                
                bhControl = (sum(semispan) * control([1 2]));
                
                % Need two to avoid control surface deflection affecting
                % non-control surface partitions
                bhControl = repmat(bhControl,1,2);
                
                [bhSort,IDorder] = sort([bh, bhControl]);
                
                controlChords = interp1(bh,chord,bhControl);
                chord = [chord, controlChords];
                chord = chord(IDorder);
                
                semispan = diff(bhSort);
                
                controlxLE = interp1(bh,xLE,bhControl);
                xLE = [xLE, controlxLE];
                xLE = xLE(IDorder);
                
                controlFoilUp = interp1(bh,FoilUp',bhControl)';
                
                FoilUp = [FoilUp, controlFoilUp];
                FoilUp = FoilUp(:,IDorder);
                
                controlFoilLo = interp1(bh,FoilLo',bhControl)';
                
                FoilLo = [FoilLo, controlFoilLo];
                FoilLo = FoilLo(:,IDorder);
                 
                % Including control spans into full config and sorting
                bh = bhSort;
                
                % Rewrite control chords to logical stating which defined
                % partitions are part of control surface
                controlSurf = bh >= bhControl(1) & bh <= bhControl(2);
                
                first = find(controlSurf,1,'first');
                last = find(controlSurf,1,'last');
                
                % Remove duplicate control surface chords from control
                % surface so that one (non-control surface) is lofted into
                % the other (begin/end control surface)
                controlSurf([first last]) = false; 
                
                [~,ID] = min(abs(control(3) - xnorm));
                
                xnorm = repmat(xnorm,1,max(IDorder));
                
                obj.Control.ChordID = ID;
                obj.Control.Surf = controlSurf;
            else
                obj.Control.ChordID = 0;
                obj.Control.Surf = zeros(nSecs,1);
            end
            
            %%
            
            z = bh * sin(di);
            
            % Calculates all xpoints with leading edge sweep offset
            [x_u,x_l] = deal(xnorm.*chord + xLE);
            
            if di < pi/4 && di > -pi/4
                y_u = bh + FoilUp * -sin(di) .* chord;
                z_u = z + FoilUp * cos(di) .* chord;
                y_l = bh + FoilLo * -sin(di) .* chord;
                z_l = z + FoilLo * cos(di) .* chord;
                
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
                obj.Name = "tail";
            else
                wing.x = [x_u, fliplr(x_l)];
                wing.y = [y_u, fliplr(y_l)];
                wing.z = [z_u, fliplr(z_l)];
                obj.Name = "wing";
            end
            
            obj.Chord = chord;
            obj.WetChord = chord;
            obj.Span = semispan;
            obj.WetSpan = semispan;
            obj.Area = Area;
            obj.WetArea = Area;
            obj.MAC = cbar;
            obj.WetMAC = cbar;
            obj.Dihedral = di;
            obj.Points = xyztopoints(wing);
            obj.Partitions = nParts;
            obj.Boolean = boolean;
            
            %%
%             plotter(obj.Points);
        end
    end
    
end

