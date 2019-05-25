function [delta,flag] = barycentrictrimodel(pt,dt,p)

x = pt(1,:);
y = pt(2,:);

Vold = pt([3 4],:);

Vcx = p(3,:);
Vcy = p(4,:);

x1 = p(1,1);
x2 = p(1,2);
x3 = p(1,3);

y1 = p(2,1);
y2 = p(2,2);
y3 = p(2,3);

W1 = ((y2 - y3).*(x - x3) + (x3 - x2).*(y - y3))./((y2 - y3).*(x1 - x3) + (x3 - x2).*(y1 - y3));
W2 = ((y3 - y1).*(x - x3) + (x1 - x3).*(y - y3))./((y2 - y3).*(x1 - x3) + (x3 - x2).*(y1 - y3));

W3 = 1 - W1 - W2;

% Reasonable extrapolation boundaries
if any([W1 W2 W3] > 1.5 | [W1 W2 W3] < -0.5)
    
    flag = true;
else
    flag = false;
end

Vx = W1.*Vcx(1) + W2.*Vcx(2) + W3.*Vcx(3);
Vy = W1.*Vcy(1) + W2.*Vcy(2) + W3.*Vcy(3);

V = [Vx; Vy];

a = (V - Vold)./dt;

delta = [V;a] * dt;