classdef wingtail
    
    properties
        Name
        Conical = false
        delta = 0.1
        Offset
        xPanels = 50
        % "Linear", "Cosine" or "HalfCosine"
        Distribution = "Cosine"
        xNorm
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
        Box
        Skin
        Spar
        Control
        TrailingEdge = false;
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
                    xNorm = X/max(X);
                case "Cosine"
                    xNorm = 0.5*(1-cos((X*pi)/max(X)));
                case "HalfCosine"
                    xNorm = 1-cos((X*(pi/2))/max(X));
            end
            
            %% Aerofoil Section & Wingbox x-Discretisation
            
            [foilUp,foilLo] = deal(zeros(numel(xNorm),nSecs));
            
            for i=1:nSecs
                
                Sec = sections{i};
                
                dim = size(Sec,1);
                half = ceil(dim/2);
                                
                upper = unique(Sec(1:half,:),'stable','rows');
                lower = unique(Sec(half:dim,:),'stable','rows');
                
                % Interpolate z given section xz coords and defined 
                % x-discretisation 
                foilUp(:,i) = interp1(upper(:,1),upper(:,2),xNorm,'pchip'); % interpolates sample data wrt defined chorwise points
                foilLo(:,i) = interp1(lower(:,1),lower(:,2),xNorm,'pchip');
            end
            
            %% Control surfaces
            if exist('control','var')
                
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
                
                controlFoilUp = interp1(bh,foilUp',bhControl)';
                
                foilUp = [foilUp, controlFoilUp];
                foilUp = foilUp(:,IDorder);
                
                controlFoilLo = interp1(bh,foilLo',bhControl)';
                
                foilLo = [foilLo, controlFoilLo];
                foilLo = foilLo(:,IDorder);
                 
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
                
                [~,ID] = min(abs(control(3) - xNorm));
                
                xNorm = repmat(xNorm,1,max(IDorder));
                
                obj.Control.ChordID = ID;
                obj.Control.Surf = controlSurf;
            else
                
                obj.Control.ChordID = 0;
                obj.Control.Surf = zeros(nSecs,1);
            end
            
            %%
            
            bh(:,:,2) = bh * sin(di);
            
            % Calculates all xpoints with leading edge sweep offset
            [wingUpper, wingLower] = deal(xNorm.*chord + xLE);
            
            if di < pi/4 && di > -pi/4
                
                % Create y, z rotation matrix, centreline chord should
                % have zero degree rotation for symmetry
                rotation = [0, 1; repmat([-sin(di) cos(di)], nParts, 1)];
                rotation = permute(rotation,[3,1,2]);
                
                wingUpper(:,:,[2 3]) = bh + foilUp .* rotation .* chord;
                wingLower(:,:,[2 3]) = bh + foilLo .* rotation .* chord;
            else
                % This doesn't work
            end
            
            if ~isequal(wingUpper(end,:,:),wingLower(end,:,:))
                
                wingUpper(end+1,:,:) = wingLower(end,:,:);
                wingLower = wingLower([1:end end],:,:);
                
                xNorm(end+1) = 1;
                obj.TrailingEdge = true;
            end
            
            if boolean
                
                points = wingLower;
                obj.Name = "tail";
            else
                
                points = [wingUpper, fliplr(wingLower)];
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
            obj.xNorm = xNorm;
            obj.Points = points;
            obj.Partitions = nParts;
            obj.Boolean = boolean;
            
            %%
%             plotter(obj.Points);
        end
    end
    
end

