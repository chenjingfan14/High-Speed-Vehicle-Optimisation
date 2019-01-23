classdef wingtail_struct
    
    properties
        Name
        Conical = false
        delta = 1
        Offset
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
            
            sparEdges = [sparBegin(1); sparEnd(end)];
            between = xnorm >= sparEdges(1) & xnorm <= sparEdges(2);
            xBox = xnorm(between);
            
            [foilUp,foilLo] = deal(zeros(numel(xnorm),nParts+1));
            [boxUp,boxLo] = deal(zeros(numel(xBox),nParts+1));
            
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
                
                boxUp(:,i) = foilUp(between,i);
                boxLo(:,i) = foilLo(between,i);
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
                
                yWingUp = bh + foilUp * -sin(di) .* chord;
                zWingUp = z + foilUp * cos(di) .* chord;
                yWingLo = bh + foilLo * -sin(di) .* chord;
                zWingLo = z + foilLo * cos(di) .* chord;
                
                yBoxUpOut = bh + boxUp * -sin(di) .* chord;
                zBoxUpOut = z + boxUp * cos(di) .* chord;
                yBoxLoOut = bh + boxLo * -sin(di) .* chord;
                zBoxLoOut = z + boxLo * cos(di) .* chord;
                
                % Interp/Extrapolate to align root with centreline
                if di~=0
                    
                    zWingUp(:,1) = zWingUp(:,1)+(-yWingUp(:,1).*((zWingUp(:,2)-zWingUp(:,1))./(yWingUp(:,2)-yWingUp(:,1))));
                    yWingUp(:,1) = 0;
                    zWingLo(:,1) = zWingLo(:,1)+(-yWingLo(:,1).*((zWingLo(:,2)-zWingLo(:,1))./(yWingLo(:,2)-yWingLo(:,1))));
                    yWingLo(:,1) = 0;
                    
                    zBoxUpOut(:,1) = zBoxUpOut(:,1)+(-yBoxUpOut(:,1).*((zBoxUpOut(:,2)-zBoxUpOut(:,1))./(yBoxUpOut(:,2)-yBoxUpOut(:,1))));
                    yBoxUpOut(:,1) = 0;
                    zBoxLoOut(:,1) = zBoxLoOut(:,1)+(-yBoxLoOut(:,1).*((zBoxLoOut(:,2)-zBoxLoOut(:,1))./(yBoxLoOut(:,2)-yBoxLoOut(:,1))));
                    yBoxLoOut(:,1) = 0;
                end
            else
                % This doesn't work
            end
            
            %% Wingbox
            
            zBoxUpIn = zBoxUpOut - (skinThick * chord);
            zBoxLoIn = zBoxLoOut + (skinThick * chord);
            
            yBoxUpIn = yBoxLoOut + (zBoxUpIn - zBoxLoOut).*((yBoxUpOut - yBoxLoOut)./(zBoxUpOut - zBoxLoOut));
            yBoxLoIn = yBoxLoOut + (zBoxLoIn - zBoxLoOut).*((yBoxUpOut - yBoxLoOut)./(zBoxUpOut - zBoxLoOut));
            
%             y = y0 + (x - x0).*((y1 - y0)./(x1 - x0));
            
            for i = numSpars:-1:1
                
                con = any(sparsFront(:,i) == xBox');
                
                ax = (sparsFront(:,i) .* chord + xLE)';
                bx = (sparsFront(:,i) .* chord + xLE)';
                
                ay = yBoxUpOut(con,:)';
                by = yBoxLoOut(con,:)';
                
                az = zBoxUpIn(con,:)';
                bz = zBoxLoIn(con,:)';
                
                [dim,~] = size(ax);
                
                spars.x = reshape([ax(:) bx(:)]', 2*dim,[])';
                spars.y = reshape([ay(:) by(:)]', 2*dim,[])';
                spars.z = reshape([az(:) bz(:)]', 2*dim,[])';
                
                spars = xyztopoints(spars);
                
                % 4 points for each spar, 3 coords per point, repeated for
                % number of sections
                sparCell(i,:) = mat2cell(spars.xyz,2,repmat(2*dim,1,nSecs));
                                
%                 xSparsFront(:,:,i) = repmat(sparsFront(i),2,1);
%                 ySparsFront(:,:,i) = [yBoxUp(con,:); yBoxLo(con,:)];
%                 zSparsFront(:,:,i) = [zBoxUpIn(con,:); zBoxLoIn(con,:)];
%                 
%                 con = sparsBack(i,:) == xBox;
%                 
%                 xSparsBack(:,:,i) = repmat(sparsBack(i),2,1);
%                 ySparsBack(:,:,i) = [yBoxUp(con,:); yBoxLo(con,:)];
%                 zSparsBack(:,:,i) = [zBoxUpIn(con,:); zBoxLoIn(con,:)];
                
                
                
            end
            
            
%             xSparsFront = xSparsFront .* chord + xLE;
%             xSparsBack = xSparsBack .* chord + xLE;
            
            xBox = xBox * chord + xLE;
            
            figure
            hold on
            axis equal
            plot3(x_u,yWingUp,zWingUp,'k')
            plot3(x_l,yWingLo,zWingLo,'k')
            
            plot3(xBox,yBoxUpOut,zBoxUpOut,'r')
            plot3(xBox,yBoxLoOut,zBoxLoOut,'r')
            plot3(xBox,yBoxUpIn,zBoxUpIn,'r')
            plot3(xBox,yBoxLoIn,zBoxLoIn,'r')
            
            plot3(xBox,yBoxUpOut,zBoxUpOut,'r.')
            plot3(xBox,yBoxLoOut,zBoxLoOut,'r.')
            plot3(xBox,yBoxUpIn,zBoxUpIn,'r.')
            plot3(xBox,yBoxLoIn,zBoxLoIn,'r.')
            
            % Stupid replication to produce entire box output
%             xSpars(end+1,:) = xSpars(1,:);
%             ySpars(end+1,:) = ySpars(1,:);
%             zSpars(end+1,:) = zSpars(1,:);
            
%             for i = 1:numSpars
%                 
%                 plot3(xSparsFront(:,i),ySparsFront(:,i),zSparsFront(:,i),'r.')
%                 plot3(xSparsBack(:,i),ySparsBack(:,i),zSparsBack(:,i),'r.')
%             end
            
            hold off
            
            %%
            if ~isequal(zWingUp(end,:),zWingLo(end,:))
                x_u(end+1,:) = x_l(end,:);
                yWingUp(end+1,:) = yWingLo(end,:);
                zWingUp(end+1,:) = zWingLo(end,:);
                x_l(end+1,:) = x_l(end,:);
                yWingLo(end+1,:) = yWingLo(end,:);
                zWingLo(end+1,:) = zWingLo(end,:);
            end
            
            if boolean
                
                wing.x = x_l;
                wing.y = yWingLo;
                wing.z = zWingLo;
                obj.Name = "tail";
            else
                
                wing.x = [x_u, fliplr(x_l)];
                wing.y = [yWingUp, fliplr(yWingLo)];
                wing.z = [zWingUp, fliplr(zWingLo)];
                obj.Name = "wing";
            end
            
            skin.Thickness = skinThick * chord;
            
            for i = nSecs:-1:1
                
                xSkinUpper(:,:,i) = repmat(xBox(:,i),1,2);
                ySkinUpper(:,:,i) = repmat(yBoxUpOut(:,i),1,2);
                zSkinUpper(:,:,i) = [zBoxUpOut(:,i), zBoxUpIn(:,i)];
                
                xSkinLower(:,:,i) = repmat(xBox(:,i),1,2);
                ySkinLower(:,:,i) = repmat(yBoxLoOut(:,i),1,2);
                zSkinLower(:,:,i) = [zBoxLoOut(:,i), zBoxLoIn(:,i)];
            end
            
            skin.Upper.x = xSkinUpper;
            skin.Upper.y = ySkinUpper;
            skin.Upper.z = zSkinUpper;
            
            skin.Lower.x = xSkinLower;
            skin.Lower.y = ySkinLower;
            skin.Lower.z = zSkinLower;
            
            spars.Thickness = sparThick * chord;
            spars.Points = sparCell;
            
%             spars.Front.x = xSparsFront;
%             spars.Front.y = ySparsFront;
%             spars.Front.z = zSparsFront;
%             
%             spars.Back.x = xSparsBack;
%             spars.Back.y = ySparsBack;
%             spars.Back.z = zSparsBack;
            
            skin.Upper = xyztopoints(skin.Upper);
            skin.Lower = xyztopoints(skin.Lower);
            
%             spars.Front = xyztopoints(spars.Front);
%             spars.Back = xyztopoints(spars.Back);
            
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
            obj.Skin = skin;
            obj.Spars = spars;
            obj.Partitions = nParts;
            obj.Boolean = boolean;
            
            %%
%             plotter(obj.Points);
        end
    end
    
end

