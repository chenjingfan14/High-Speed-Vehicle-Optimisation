function within = iswithintriangle(P,c)

% Barycentric

% x = pt(1,:);
% y = pt(2,:);
% 
% x1 = p(1,1);
% x2 = p(1,2);
% x3 = p(1,3);
% 
% y1 = p(2,1);
% y2 = p(2,2);
% y3 = p(2,3);
% 
% W1 = ((y2 - y3).*(x - x3) + (x3 - x2).*(y - y3))./((y2 - y3).*(x1 - x3) + (x3 - x2).*(y1 - y3));
% W2 = ((y3 - y1).*(x - x3) + (x1 - x3).*(y - y3))./((y2 - y3).*(x1 - x3) + (x3 - x2).*(y1 - y3));
% 
% W3 = 1 - W1 - W2;
% 
% W = [W1 W2 W3];
% 
% within = all(W > 0 & W < 1);

% W. Heidrich, Journal of Graphics, GPU, and Game Tools,Volume 10, Issue 3, 2005

P1 = c(:,1);
P2 = c(:,2);
P3 = c(:,3);

u = P2 - P1;
v = P3 - P1;
w = P - P1;

n = crossmat(u,v);

gamma = (crossmat(u,w).*n)/(n.^2);
beta = (crossmat(w,v).*n)/(n.^2);

alpha = 1 - gamma - beta;

abg = [alpha beta gamma];

within = all(abg >= 0 & abg <= 1);

% https://math.stackexchange.com/questions/4322/check-whether-a-point-is-within-a-3d-triangle

% A = crossmat(P2 - P1, P3 - P1);
% C = crossmat(P3 - P, P1 - P);
% D = crossmat(P1 - P, P2 - P);
% 
% k1 = C/A;
% k2 = D/A;
% 
% within = all([k1 k2] >= -0.1);

