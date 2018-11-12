classdef arbitraryfuse < body
    %CIRCLEFUSELAGE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        delta = 0.1
        Height
        Width
        Geometry
        Length
        Area
    end
    
    methods
        function obj=arbitraryfuse(input,offset)
            
            % [Upper L, Top Smaj]
            yTop = [input(1),input(2)];
            
            yTot = sum(yTop);
            % [Bot Smaj, Bot L]
            yBot=[input(3)*yTot,(1-input(3))*yTot];
            % [Top Smin, Side L, Bot Smin]
            z=[input(4),input(5),input(6)];
            
            aftLength = input(7);
            
            delta = obj.delta; % meters
            
            width = yTot;
            height = sum(z);
            geo = [yTop,z(2),yBot];
            smin = z([1,3]);
            
            [dim1,dim2] = size(geo);
            [A,numseg] = deal(zeros(dim1,dim2));
            
            for i=1:dim2
                
                igeo = geo(i);
                
                switch i
                    case {1,5}
                        numseg(i) = ceil(igeo/delta);
                        A(i) = igeo*height/2;
                    case 3
                        numseg(i) = ceil(igeo/delta);
                        A(i) = igeo*width;
                    case {2,4}
                        a = igeo;
                        b = smin(i/2);
                        
                        norm = sqrt(a^2 + b^2);
                        if norm < delta*5
                            idelta = norm/5;
                        else
                            idelta = delta;
                        end
                        
                        h = ((a-b)^2)/((a+b)^2);
                        p = 0.25*pi*(a+b)*(1+(3*h/(10+sqrt(4-3*h))));
                        numseg(i) = ceil(p/idelta);
                        A(i) = 0.25*pi*a*b;
                end
            end
            
            dim3 = max(numseg)+1;
            y = NaN(dim2,dim3);
            z = NaN(dim2,dim3);
            
            for i=1:dim2
                
                igeo = geo(i);
                iseg = numseg(i);
                
                % If section geometry is zero, or semiminor axis for radius
                % sections is zero, continue
                if igeo==0 || (ismember(i,[2,4]) && smin(i/2) == 0)
                    continue
                else
                    
                    switch i
                        case 1 % Top straight
                            
                            fill = 0:igeo/iseg:igeo;
                            [~,dim] = size(fill);
                            cols = 1:dim;
                            
                            y(i,cols) = fill;
                            z(i,cols) = 0;
                            
                        case 2 % Top rad
                            
                            a = geo(i);
                            b = smin(1);
                            multi = 0:iseg;
                            cols = multi+1;
                            
                            y(i,cols) = a*sin(multi*(pi/2)/iseg);
                            z(i,cols) = b*cos(multi*(pi/2)/iseg);
                            
                        case 3 % Side straight
                            
                            fill = igeo:-igeo/iseg:0;
                            [~,dim] = size(fill);
                            cols = 1:dim;
                            
                            z(i,cols) = fill;
                            y(i,cols) = 0;
                            
                        case 4 % Bottom rad
                            
                            a = geo(i);
                            b = smin(2);
                            multi = 0:iseg;
                            cols = multi+1;
                            
                            y(i,cols) = a*cos(multi*(pi/2)/iseg);
                            z(i,cols) = -b*sin(multi*(pi/2)/iseg);
                            
                        case 5 % Bottom straight

                            fill = igeo:-igeo/iseg:0;
                            [~,dim] = size(fill);
                            cols = 1:dim;
                            
                            y(i,cols) = fill;
                            z(i,cols) = 0;
                    end
                end
            end
            
            con = ~isnan(y);
            keep = any(con,2);
            
            con = con(keep,:);
            y = y(keep,:);
            z = z(keep,:);
            
            s.x = [0;aftLength] + offset;
            
            rows = sum(keep);
            tots = sum(con,2);
            tot = sum(tots)-rows+1;
            
            [yint,zint] = deal(zeros(1,tot));
            last = 1;
            
            for i=1:rows
                
                coni = con(i,:);
                idim = tots(i);
                cols = (0:idim-1)+last;
                
                if i==1
                    yint(cols) = y(i,coni);
                    zint(cols) = z(i,coni);
                else
                    diffy = yint(last)-y(i,1);
                    diffz = zint(last)-z(i,1);
                
                    yi = y(i,coni) + diffy;
                    zi = z(i,coni) + diffz;

                    yint(cols) = yi;
                    zint(cols) = zi;
                end
                
                last = cols(end);
            end
            
            zint = zint-((zint(1)+zint(end))/2);
            yint(end) = 0;
            
            s.y = yint;
            s.z = zint;
            
            obj.Area = sum(A);
            obj.Height = height;
            obj.Width = width*2;
            obj.Length = aftLength;
            obj.Geometry = geo;
            obj.Points = xyztopoints(s);
        end
    end
end

