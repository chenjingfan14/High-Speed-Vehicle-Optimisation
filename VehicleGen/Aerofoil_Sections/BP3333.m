% function foilData = Bezier3(controlPointsArray,nParts,nBez,nPop)

%Bezier Curves
%Script to plot an aerofoil section using Bezier Curves

nPop = 1;
nParts = 1;

varMin = [0.05 0.01 0.1 0.2 0 -0.2 0 0 0.05 -0.04 0 0.15 0.05 -0.5 0 0 0.001];
varMax = [0.1 0.1 0.3 0.5 0.2 0 0.9 0.01 0.1 -0.001 0.7 0.4 0.15 0 0.9 0.001 0.3];

nVar = length(varMin);
nBez = 7;

varMinMat = repmat(varMin,nPop,1);
varMaxMat = repmat(varMax,nPop,1);

varSize = [nPop nVar];

count = 0;
success = false;
while ~success

    count = count + 1;
    
variables = unifrnd(varMinMat,varMaxMat,varSize);

gle = variables(:,1);
b0 = variables(:,2);
b2 = variables(:,3);
xc = variables(:,4);
yc = variables(:,5);
kc = variables(:,6);
b17 = variables(:,7);
zte = variables(:,8);
ate = variables(:,9);
rle = variables(:,10);
b8 = variables(:,11);
xt = variables(:,12);
yt = variables(:,13);
kt = variables(:,14);
b15 = variables(:,15);
dzte = variables(:,16);
Bte = variables(:,17);

[~,total] = size(variables);
array = 1:total;

% nSecs = nParts + 1;
nSecs = 1;

totPerSec = total/nSecs;

mat = reshape(array,totPerSec,nSecs);

foilData = cell(nPop,nSecs);

u = (0:1/99:1)';

%% BP3333 Paramters

a = 27/4 * kt.^2;
b = -27 * kt.^2 * xt;
c = (9 * kt * yt) + (81/2 * kt.^2 * xt.^2);
d = (2 * rle) - (18 * kt * xt * yt) - (27 * kt.^2 * xt.^3);
e = (3 * yt.^2) + (9 * kt * xt.^2 * yt) + (27/4 * kt.^2 * xt.^4);

r = real(roots([a b c d e]));

id = max(0, xt - ((-2 * yt)/(3 * kt)).^0.5) < r & r < xt;

if any(id)
    
    b9 = min(r(id));
else
    continue
end

cgle = cot(gle);
cate = cot(ate);

a = 1/(3*kc*(cgle + cate).^2);
b = 16 + 3*kc*(cgle + cate) * (1 + zte * cate);
c = 4 * (16 + 6*kc*(cgle + cate) * (1 - yc*(cgle + cate) + zte * cate)).^0.5;

b1 = [a * (b + c), a * (b - c)];

a = (16 + 3*kc*(cgle + cate) * (1 + zte * cate))/(3*kc*(cgle + cate));
b = 4 * (16 + 6*kc*(cgle + cate) * (1 - yc*(cgle + cate) + zte * cate)).^0.5;

b1 = [a + b, a - b];

id = 0 < b1 & b1 < yc;

if any(id) && isreal(b)
    
    [~,id] = min(yc - b1(id));
    b1 = b1(id);
else
    continue
end
% Leading thickness curve
% xlt = [0 0 b9 xt];

xlt0 = 0;
xlt1 = 0;
xlt2 = b9;
xlt3 = xt;

ylt0 = 0;
ylt1 = 1.5 * kt * (xt - b9).^2 + yt;
ylt2 = yt;
ylt3 = yt;

% ylt = [0 ylt1 yt yt];

% Leading camber curve
xlc0 = 0;
xlc1 = b1 * cgle;
xlc2 = xc - (2*(b1 - yc)/(3*kc)).^0.5;
xlc3 = xc;

ylc0 = 0;
ylc1 = b1;
ylc2 = yc;
ylc3 = yc;

xlc = [0 xlc1 xlc2 xc];
ylc = [0 b1 yc yc];

% Trailing thickness curve
xtt0 = xt;
xtt1 = 2*xt - b9;
xtt2 = 1 + (dzte - (1.5*kt*(xt - b9).^2 + yt)) * cot(Bte);
xtt3 = 1;

xtt = [xt xtt1 xtt2 1];

ytt0 = yt;
ytt1 = yt;
ytt2 = 1.5*kt*(xt - b9).^2 + yt;
ytt3 = dzte;

ytt = [yt yt ytt2 dzte];

% Trailing camber curve
xtc0 = xc;
xtc1 = xc + (2*(b1 - yc)/(3 * kc)).^0.5;
xtc2 = 1 + (zte - b1) * cate;
xtc3 = 1;

ytc0 = yc;
ytc1 = yc;
ytc2 = b1;
ytc3 = zte;

xtc = [xc xtc1 xtc2 1];
ytc = [yc yc b1 zte];

Pxlt = xlt0*(1 - u).^3 + 3*xlt1*u.*(1 - u).^2 + 3*xlt2*u.^2.*(1 - u) + xlt3*u.^3;
Pxlc = xlc0*(1 - u).^3 + 3*xlc1*u.*(1 - u).^2 + 3*xlc2*u.^2.*(1 - u) + xlc3*u.^3;
Pxtt = xtt0*(1 - u).^3 + 3*xtt1*u.*(1 - u).^2 + 3*xtt2*u.^2.*(1 - u) + xtt3*u.^3;
Pxtc = xtc0*(1 - u).^3 + 3*xtc1*u.*(1 - u).^2 + 3*xtc2*u.^2.*(1 - u) + xtc3*u.^3;

Pylt = ylt0*(1 - u).^3 + 3*ylt1*u.*(1 - u).^2 + 3*ylt2*u.^2.*(1 - u) + ylt3*u.^3;
Pylc = ylc0*(1 - u).^3 + 3*ylc1*u.*(1 - u).^2 + 3*ylc2*u.^2.*(1 - u) + ylc3*u.^3;
Pytt = ytt0*(1 - u).^3 + 3*ytt1*u.*(1 - u).^2 + 3*ytt2*u.^2.*(1 - u) + ytt3*u.^3;
Pytc = ytc0*(1 - u).^3 + 3*ytc1*u.*(1 - u).^2 + 3*ytc2*u.^2.*(1 - u) + ytc3*u.^3;

success = true;

figure
hold on
axis equal
plot(Pxlt,Pylt,'r');
plot(Pxlc,Pylc,'b');
plot(Pxlt,Pylt,'r');
plot(Pxtc,Pytc,'b');
xlabel('x/c');
ylabel('y/c');
title('Aerofoil Section');

end

count

for i = 1:nSecs
    
    split = mat(:,i);
    
    for j = 1:nPop
        %set i for required number of aerofoil configurations
        % upper x coordinates of control points
        
        x = [0 variables(j,split([1 3 5 6 7])) 1];
        y = [0 variables(j,split([2 4 4 4 8])) 0];
        xt = [0 0 variables(j,split([10 12 13 14])) 1];
        t = [0 variables(j,split([9 11 11 11 15])) 0];
        %check number of x co-ordinates is equal to number of y co-ordinates
        
        b = f.*((1-u).^(n-k)).*(u.^k).*x;
        c = f.*((1-u).^(n-k)).*(u.^k).*y;
        
        Px = sum(b,2);
        Pz = sum(c,2);
        
        b = f.*((1-u).^(n-k)).*(u.^k).*xt;
        c = f.*((1-u).^(n-k)).*(u.^k).*t;
        
        Pxt = sum(b,2);
        Pt = sum(c,2);
                 
        xdiff = diff(Px);
        zdiff = diff(Pz);
        
        xdiff(end+1) = -(Px(end-1) - Px(end));
        zdiff(end+1) = Pz(end-1) - Pz(end); 
        
        theta = atan(abs(zdiff)./xdiff);
        
        theta(isnan(theta)) = 0;
        
        con = zdiff < 0;
        theta(con) = -theta(con);
        
        xu = Px + Pt .* cos(theta);% - t .* sin(theta);
        zu = Pz + Pt .* sin(theta);% + t .* sin(theta);
        
        xl = Px - Pt .* cos(theta);% - t .* sin(theta);
        zl = Pz - Pt .* sin(theta);
        
        % plot aerofoil
        figure
        hold on
%         axis equal
        plot(Px,Pz,'r');
        plot(x,y,'rx');
        plot(xu,zu,'b');
        plot(xl,zl,'b');
        xlabel('x/c');
        ylabel('y/c');
        title('Aerofoil Section');
%         axis([0 1 -0.1 0.1]);
        legend({'Camber','Control Points'},'Location','northeast');
        hold off
        
%         foilData{j,i} = [Px, Pz; flipud([Pxl, Pzl])];
        
    end
end