function [intersection,reason] = planeintersection(n,V0,P0,P1,xBody,yzBody,j)

intersection = [];
reason = [];

u = P1 - P0;
w = P0 - V0;
D = n(:,1).*u(:,1) + n(:,2).*u(:,2) + n(:,3).*u(:,3);
N = -(n(:,1).*w(:,1) + n(:,2).*w(:,2) + n(:,3).*w(:,3));

% Compute the Isection parameter
sI = N ./ D;
I = P0 + sI*u;

yBody = yzBody(:,1);
zBody = yzBody(:,2);

conx = I(:,1) >= xBody(1) & I(:,1) <= xBody(2) | I(:,1) >= xBody(2) & I(:,1) <= xBody(1);
cony = I(:,2) >= yBody(j) & I(:,2) <= yBody(j+1) | I(:,2) >= yBody(j+1) & I(:,2) <= yBody(j);
conz = I(:,3) >= zBody(j) & I(:,3) <= zBody(j+1) | I(:,3) >= zBody(j+1) & I(:,3) <= zBody(j);

inter = I(conx & cony & conz,:);

% Sometimes needs to be rounded, so if first attempt returns empty, round
if isempty(inter)
    
    I = round(I,10);
    
    conx = I(:,1) >= xBody(1) & I(:,1) <= xBody(2) | I(:,1) >= xBody(2) & I(:,1) <= xBody(1);
    cony = I(:,2) >= yBody(j) & I(:,2) <= yBody(j+1) | I(:,2) >= yBody(j+1) & I(:,2) <= yBody(j);
    conz = I(:,3) >= zBody(j) & I(:,3) <= zBody(j+1) | I(:,3) >= zBody(j+1) & I(:,3) <= zBody(j);

    inter = I(conx & cony & conz,:);
end

if isempty(inter)

    % Initial attempt at finding problem
    if any(I(cony & conz,1) <= xBody(1))
        % Intersection before body
        reason = 1;
    elseif any(I(cony & conz,1) >= xBody(2))
        % Intersection beyond body
        reason = 2;
    elseif any(I(conx,3) >= zBody(j(1)))
        % Intersection above body
        reason = 3;
    elseif any(I(conx,3) <= zBody(j(end)))
        % Intersection below body
        reason = 4;
        
    elseif any((P1(2).^2 + (P1(3) - P0(3)).^2).^0.5 < (yBody(j).^2 + zBody(j).^2).^0.5)
        %% CHECK THIS DOES THE RIGHT THING
        % Check that end of first partition isn't inside body
        % Increase span
        reason = 5;
    
    % Second (less refined) attempt
    elseif P0(3) <= zBody(j(1))
        reason = 3;
    elseif P0(3) >= zBody(j(end))
        reason = 4;
    elseif P0(1) <= xBody(1)
        reason = 1;
    elseif P0(1) >= xBody(2)
        reason = 2;
    else
        reason = 5;
    end

else
    % Ensures intersection point is between P0 and P1, outwith will return
    % a infeasible design
    between = (inter(:,1) >= P0(1) & inter(:,1) <= P1(1)) | (inter(:,1) >= P1(1) & inter(:,1) <= P0(1));
    
    if between
        % Ensure only one intersection is returned
        intersection = inter(between(1),:);
    else
        % Increase semispan
        reason = 5;
    end
end