classdef nose < body
    %NOSE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name = "nose"
        Radius
        Length
        zOff
        xPanels
        Theta
        delta
        Rotation
    end
    
    methods
        function obj=nose(input,dim_f,NFLength)
            
            r = input(1);
            ln = input(2);
            zOff = input(3);
            delta = pi/(dim_f-1);
            xPanels = 10;
            phi = 0:delta:pi;
            
            if all([r,ln])
                
                rot = abs(atan(zOff/NFLength));
                rotArm = (abs(zOff)^2 + NFLength^2)^0.5;
                thetaMax = acos(-(ln/r)+1);
                theta = (0:thetaMax/xPanels:thetaMax)';
                
                x = (r*(1-cos(theta))*ones(1,length(phi))) + rotArm;
                y = r*sin(theta).*sin(phi);
                z = r*cos(phi).*sin(theta);
                
            else
                
                rot = 0;
                rotArm = 0;
                theta = 0;
                dim = size(phi);
                [x,y,z] = deal(zeros(dim));
      
            end
            
            x = (x*cos(rot) - z*sin(rot)) - rotArm;
            a.x = x - min(x(:));
            a.y = y;
            a.z = z*cos(rot)+(x-rotArm)*sin(rot);
            obj.Points = xyztopoints(a);
            
            obj.Radius = r;
            obj.Length = ln;
            obj.zOff = zOff;
            obj.xPanels = xPanels;
            obj.Theta = theta;
            obj.delta=delta;
            obj.Rotation=rot;
        end
        
        function plot(obj)
            figure
            hold on
            axis equal
            xlabel('x')
            ylabel('y')
            zlabel('z')
            [x,y]=size(obj.Points);
            a=obj.Points;
            for i=1:x
                for j=1:y
                    plot3(a{i,j}(1),a{i,j}(2),a{i,j}(3),'k*')
                end
            end
        end
    end
end