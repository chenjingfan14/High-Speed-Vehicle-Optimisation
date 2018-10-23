clear all

quarcirc = 0:pi/5:pi/2;
rad = 1;

xfqc = (-rad*sin(quarcirc))+rad;
yfqc = zeros(size(xfqc));
zfqc = rad*cos(quarcirc);

xbqc_off = 7;
ybqc_off = 5;

xbqc = ones(size(xfqc))*xbqc_off;
ybqc = xfqc - ybqc_off;
zbqc = zfqc;

xPanels = 20;

for i = 1:length(xfqc)
    xbegin = xfqc(i);
    xfinish = xbqc(i);
    x(i,:) = xbegin:(xfinish-xbegin)/xPanels:xfinish;
    y(i,:)=yfqc(i)+(x(i,:)-xbegin).*((ybqc(i)-yfqc(i))./(xfinish-xbegin));
    z(i,:)=zfqc(i)+(x(i,:)-xbegin).*((zbqc(i)-zfqc(i))./(xfinish-xbegin));
end

figure
hold on
axis equal
plot3(x,y,z,'k*')
xlabel('x')
ylabel('y')
zlabel('z')