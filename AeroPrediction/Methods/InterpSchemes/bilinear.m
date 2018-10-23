function dot = bilinear(p,dt,pV,V1,V2,V3,V4)

xy = p(1:2,:);
Vold = p(3:4,:);

% Implement properly for arbitrary rectangle
% u = x./max(x(:)) - x(1);
% v = y./max(y(:)) - y(1);

diffx = xy - pV;

u = diffx(1,:);
v = diffx(2,:);

v0 = (1-u).*V1 + u.*V2;
v1 = (1-u).*V4 + u.*V3;

V = (1-v).*v0 + v.*v1;

a = (V - Vold)./dt;

dot = [V;a];