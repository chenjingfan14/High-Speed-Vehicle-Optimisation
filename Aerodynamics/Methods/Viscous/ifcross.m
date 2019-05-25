function [crossCoords, cross, p2, Vp2] = ifcross(p1,p2,Vp2,c1,c2,V1,V2,crossCoords)

int = lineintersection3D(p1,p2,c1,c2);

% Ensure intersection is with line segment
within = all((int >= p1 & int <= p2 | int <= p1 & int >= p2) &...
            (int >= c1 & int <= c2 | int <= c1 & int >= c2));

if within
    
    % Which panel line is being crossed
    test = [c1 c2];
    
    % Particle cannot cross back over line it has just crossed, must be
    % different otherwise no intersection
    if isequal(test, crossCoords) || isequal(test(:,[2 1]), crossCoords)
        
        cross = false;
    else
        crossCoords = test;
        p2 = int;
        Vp2 = V1 + (p2 - c1).*((V2 - V1)./(c2 - c1));
        
        % NaN values come from c1, c2 x, y or z values are equal. In such
        % cases, an average value is taken between V1, V2
        con = isnan(Vp2);
        
        if any(con)
            
            Vp2(con) = (V1(con) + V2(con))/2;
        end
        
        cross = true;
    end
else
    cross = false;
end