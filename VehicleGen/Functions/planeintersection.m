function [intersection,reason] = planeintersection(n,V0,P0,P1,xBody,yBody,zBody,j)

intersection = [];
reason = [];

dim = length(V0);

X = 1:3:dim;
Y = X + 1;
Z = Y + 1;

matV0 = [V0(X),V0(Y),V0(Z)];

u = P1 - P0;
w = P0 - matV0;
D = n(:,1).*u(:,1) + n(:,2).*u(:,2) + n(:,3).*u(:,3);
N = -(n(:,1).*w(:,1) + n(:,2).*w(:,2) + n(:,3).*w(:,3));

% Compute the Isection parameter
sI = N ./ D;
I = P0 + sI*u;

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

    if any([P0(1);I(:,1)] <= xBody(1))
        % Intersection before body
        reason = 1;
    elseif any([P0(1);I(:,1)] >= xBody(2))
        % Intersection beyond body
        reason = 2;
    elseif any([P0(3);I(:,3)] >= zBody(j(1)))
        % Intersection above body
        reason = 3;
    elseif any([P0(3);I(:,3)] <= zBody(j(end)))
        % Intersection below body
        reason = 4;
    else
        % Increase span
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