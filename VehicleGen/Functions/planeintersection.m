function intersection = planeintersection(n,V0,P0,P1,xBody,yBody,zBody,j)

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

I = I(conx & cony & conz,:);

if isempty(I)
    intersection = [];
else
    
    % Ensures intersection point is between P0 and P1, outwith will return
    % a infeasible design
    between = (I(1) >= P0(1) & I(1) <= P1(1)) | (I(1) >= P1(1) & I(1) <= P0(1));

    % Ensure only one intersection is returned
    intersection = I(between(1),:);
end