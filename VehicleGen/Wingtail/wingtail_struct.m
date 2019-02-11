classdef wingtail_struct
    
    properties
        Name
        Conical = false
        delta = 1
        Offset
        % CHANGE: setting as +1
        xPanels = 5
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
        Box
        Skin
        Spars
        Control
        Boolean
        
    end
    
    methods
        function obj=wingtail_struct(dihedral,semispan,chord,sweep,sections,wingbox,control)
            
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
            
            %% Aerofoil Section & Wingbox x-Discretisation
            
            sparMid = wingbox.Location;
            skinThick = wingbox.SkinThick;% * chord;
            sparThick = wingbox.SparThick;% * chord;
            
            sparBegin = sparMid - sparThick/2;
            sparEnd = sparMid + sparThick/2;
            
            numSpars = numel(sparBegin);
            
            xnorm = sort([xnorm; sparBegin; sparEnd]);
            
            for i = numSpars:-1:1
            
                con = xnorm >= sparBegin(i) & xnorm <= sparEnd(i);
                
                first = find(con,1,'first');
                last = find(con,1,'last');
                
%                 sparsFront(i,:) = xnorm(first);
%                 sparsBack(i,:) = xnorm(last);
                
                sparsFront(:,i) = xnorm([first last]);

                if sum(con) > 2
                    
                    delete = first + 1:last - 1;
                    
                    xnorm(delete) = [];
                end
            end
            
            sparsFront = sparsFront(:);
            
            sparEdges = [sparBegin(1); sparEnd(end)];
            between = xnorm >= sparEdges(1) & xnorm <= sparEdges(2);
            xBox = xnorm(between);
            
            [foilUp,foilLo] = deal(zeros(numel(xnorm),nSecs));
            
            for i=1:nSecs
                
                Sec = sections{i};
                
                dim = size(Sec,1);
                half = ceil(dim/2);
                                
                upper = unique(Sec(1:half,:),'stable','rows');
                lower = unique(Sec(half:dim,:),'stable','rows');
                
                % Interpolate z given section xz coords and defined 
                % x-discretisation 
                foilUp(:,i) = interp1(upper(:,1),upper(:,2),xnorm,'pchip'); % interpolates sample data wrt defined chorwise points
                foilLo(:,i) = interp1(lower(:,1),lower(:,2),xnorm,'pchip');
            end
               
            boxUp = foilUp(between,:);
            boxLo = foilLo(between,:);
            
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
                
                [~,ID] = min(abs(control(3) - xnorm));
                
                xnorm = repmat(xnorm,1,max(IDorder));
                
                obj.Control.ChordID = ID;
                obj.Control.Surf = controlSurf;
            else
                
                obj.Control.ChordID = 0;
                obj.Control.Surf = zeros(nSecs,1);
            end
            
            %%
            
            bh(:,:,2) = bh * sin(di);
            
            % Calculates all xpoints with leading edge sweep offset
            [wingUpper, wingLower] = deal(xnorm.*chord + xLE);
            
            if di < pi/4 && di > -pi/4
                
                % Create y, z rotation matrix, centreline chord should
                % have zero degree rotation for symmetry
                rotation = [0, 1; repmat([-sin(di) cos(di)], nParts, 1)];
                rotation = permute(rotation,[3,1,2]);
                
                wingUpper(:,:,[2 3]) = bh + foilUp .* rotation .* chord;
                wingLower(:,:,[2 3]) = bh + foilLo .* rotation .* chord;
                
                yzBoxUpperOut = bh + boxUp .* rotation .*chord;
                yzBoxLowerOut = bh + boxLo .* rotation .*chord;
                
                yzBoxUpperIn = yzBoxUpperOut - (skinThick * rotation .* chord);
                yzBoxLowerIn = yzBoxLowerOut + (skinThick * rotation .* chord);
            else
                % This doesn't work
            end
            
            [boxUpperOut, boxLowerOut, boxUpperIn, boxLowerIn] = deal(xBox * chord + xLE); 
            boxUpperOut(:,:,[2 3]) = yzBoxUpperOut;
            boxLowerOut(:,:,[2 3]) = yzBoxLowerOut;
            boxUpperIn(:,:,[2 3]) = yzBoxUpperIn;
            boxLowerIn(:,:,[2 3]) = yzBoxLowerIn;
            
            %% Wingbox
                
            con = any(sparsFront == xBox');
            
            a = boxUpperIn(con,:,:);
            b = boxLowerIn(con,:,:);
            
            [row,~] = size(a);
            
            j = 2*row;
            for i = 3:-1:1
                
                sparPoints(:,:,i) = reshape([a(:,:,i) b(:,:,i)]',[],j)';
            end
            
            for i = numSpars:-1:1
                
                sparCell{i,:} = sparPoints(j-3:j,:,:);
                j = j - 4;
            end
            
            %% Proof plotting
%             figure
%             hold on
%             axis equal
%             plot3(wingUpper(:,:,1),wingUpper(:,:,2),wingUpper(:,:,3),'k')
%             plot3(wingLower(:,:,1),wingLower(:,:,2),wingLower(:,:,3),'k')
%             
%             plot3(boxUpperOut(:,:,1),boxUpperOut(:,:,2),boxUpperOut(:,:,3),'r')
%             plot3(boxLowerOut(:,:,1),boxLowerOut(:,:,2),boxLowerOut(:,:,3),'r')
%             
%             plot3(boxUpperIn(:,:,1),boxUpperIn(:,:,2),boxUpperIn(:,:,3),'r')
%             plot3(boxLowerIn(:,:,1),boxLowerIn(:,:,2),boxLowerIn(:,:,3),'r')
% 
%             for i = 1:numSpars
%                 
%                 plot3(sparCell{i}(:,:,1),sparCell{i}(:,:,2),sparCell{i}(:,:,3),'r.')
%                 plot3(sparCell{i}(:,:,1),sparCell{i}(:,:,2),sparCell{i}(:,:,3),'r.')
%             end
%             
%             hold off
            
            %%
            if ~isequal(wingUpper(end,:,:),wingLower(end,:,:))
                
                wingUpper(end+1,:,:) = wingLower(end,:,:);
                wingLower = wingLower([1:end end],:,:);
                
                between(end+1) = false;
                
            end
            
            if boolean
                
                points = wingLower;
                obj.Name = "tail";
            else
                
                points = [wingUpper, fliplr(wingLower)];
                obj.Name = "wing";
            end
            
            skin.Thickness = skinThick * chord;
            
            for i = nSecs:-1:1
                
                skin.Upper{i} = [boxUpperOut(:,i,:), boxUpperIn(:,i,:)];
                skin.Lower{i} = [boxLowerOut(:,i,:), boxLowerIn(:,i,:)];
                skin.UpperMean{i} = squeeze((boxUpperOut(:,i,:) + boxUpperIn(:,i,:))/2);
                skin.LowerMean{i} = squeeze((boxLowerOut(:,i,:) + boxLowerIn(:,i,:))/2);
            end
            
            spars.Thickness = sparThick * chord;
            spars.Points = sparCell;
            
            obj.Chord = chord;
            obj.WetChord = chord;
            obj.Span = semispan;
            obj.WetSpan = semispan;
            obj.Area = Area;
            obj.WetArea = Area;
            obj.MAC = cbar;
            obj.WetMAC = cbar;
            obj.Dihedral = di;
            obj.Points = points;
            obj.Box = between;
            obj.Skin = skin;
            obj.Spars = spars;
            obj.Partitions = nParts;
            obj.Boolean = boolean;
            
            %%
%             plotter(obj.Points);
        end
    end
    
end

