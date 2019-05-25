function int = lineintersection3D(P1,P2,P3,P4)

% https://math.stackexchange.com/questions/270767/find-intersection-of-two-3d-lines/270793
% P1 = [10 10 6];
% P2 = [5 5 4];
% P3 = [10 10 3];
% P4 = [5 5 5];

% P1 = [12 15 4];
% P2 = [6 8 4];
% P3 = [12 15 6];
% P4 = [6 8 2];

a = P1 - P2;
b = P3 - P4;

c = cross(a,b);

% const1 = dot(P2,c);
% const2 = dot(P4,c);
% 
% if const1 == const2
%     
%     intersect = true;
% end

d = cross(a,c);

const1 = dot(P2,d);

t = (const1 - dot(d,P4))/dot(d,b); 

int = P4 + b * t;