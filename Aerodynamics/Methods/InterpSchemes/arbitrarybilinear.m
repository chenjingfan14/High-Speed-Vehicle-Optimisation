%% DOESN'T WORK

clear all
close all

%% Creating uniform/randomly perturbed panels
% Replace with size(points.x)
x = ones(4); % x,y points (x-1,y-1 panels)

% Creating x,y point matrices from 0 to [x y]-1
% Number of points
[dim1,dim2] = size(x);

% Number of panels
dim3 = dim1-1;
dim4 = dim2-1;

points(:,:,1) = repmat((0:dim1-1)',1,dim2);
points(:,:,2) = repmat((0:dim2-1),dim1,1);
points(:,:,3) = zeros(dim1,dim2);

%% Random perturbation to make uniform grid chaotic, only inner points
randx = rand(dim3-1,dim4-1)*0.5 - 0.25;
randy = rand(dim3-1,dim4-1)*0.5 - 0.25;

points(2:end-1,2:end-1,1) = points(2:end-1,2:end-1,1) + randx;
points(2:end-1,2:end-1,2) = points(2:end-1,2:end-1,2) + randy;

plotter(points)

%% Start of interp

x1 = points(1:end-1,2:end,1);
x2 = points(2:end,2:end,1);
x3 = points(1:end-1,1:end-1,1);
x4 = points(2:end,1:end-1,1);

y1 = points(1:end-1,2:end,2);
y2 = points(2:end,2:end,2);
y3 = points(1:end-1,1:end-1,2);
y4 = points(2:end,1:end-1,2);

Px = (x1+x2+x3+x4)/4;
Py = (y1+y2+y3+y4)/4;

figure(1)
hold on
plot(Px,Py,'r*')
hold off

P1 = [50 50];
P2 = [50 0];
P3 = [50 50];
P4 = [50 0];

x31 = abs(x3 - x1);
y31 = abs(y3 - y1);
x21 = abs(x2 - x1);
y21 = abs(y2 - y1);
x42 = abs(x4 - x2);
y42 = abs(y4 - y2);
x43 = abs(x4 - x3);
y43 = abs(y4 - y3);

% Horizontal/Chordwise panels are parallel
con1 = (y4 - y3).*(x2 - x1) ~= (y2 - y1).*(x4 - x3);
con2 = (y4 - y2).*(x3 - x1) == (y3 - y1).*(x4 - x2) & ~con1;
con3 = ~con1 & ~con2;

t = zeros(dim3,dim4);
s = zeros(dim3,dim4);

%% Start conditioning here
% If arbitrary or chordwise lines parallel
if any(con1)
    A = x31.*y42 - y31.*x42;
    B = Py.*(x42 - x31) - Px.*(y42 - y31) + x31.*y2 - y31.*x2 + x1.*y42 - y1.*x42;
    C = Py.*x21 - Px.*y21 + x1.*y2 - x2.*y1;

    t1 = (-B + ((B.^2) - (4*A.*C)).^0.5)./(2*C);
    t2 = (-B - ((B.^2) - (4*A.*C)).^0.5)./(2*C);

    one = t1 >= 0 & t1 <= 1;
    two = t2 >= 0 & t2 <= 1;
    con1 = one | two;
    con2 = ~con1;

    t(con1) = real([t1(one & con1);t2(two & con1)]);

    Ay = y1(con1) + y31(con1).*t(con1);
    By = y2(con1) + y42(con1).*t(con1);
    Ax = x1(con1) + x31(con1).*t(con1);
    Bx = x2(con1) + x42(con1).*t(con1);
    
    s(con1) = (Py(con1) - Ay)./(By - Ay);
    
    figure(1)
    hold on
    plot(Ax,Ay,'b*')
    hold off
end

%%
if any(con2)
    A = x21.*y43 - y21.*y43;
    B = Py.*(x43 - x21) - Px.*(y43 - y21) + x1.*y43 - y1.*x43 + x21.*y3 - y21.*x3;
    C = Py.*x31 - Px.*y31 + x1.*y3 - x3.*y1;

    s1 = (-B + ((B.^2) - (4.*A.*C)).^0.5)./(2.*C);
    s2 = (-B - ((B.^2) - (4.*A.*C)).^0.5)./(2.*C);

    one = s1 >= 0 & s1 <= 1;
    two = s2 >= 0 & s2 <= 1;

    s(con2) = real([s1(one & con2);s2(two & con2)]);
    
    Ay = y1(con2) + y31(con2).*t(con2);
    By = y2(con2) + y42(con2).*t(con2);
    Ax = x1(con2) + x31(con2).*t(con2);
    Bx = x2(con2) + x42(con2).*t(con2);
    
    Cx = x1(con2) + x21(con2).*s(con2);
    Dx = x3(con2) + x43(con2).*s(con2);
    Cy = y1(con2) + y21(con2).*s(con2);
    Dy = y3(con2) + y43(con2).*s(con2);

    t(con2) = (Py(con2) - Cy)./(Dy - Cy);
    
    figure(1)
    hold on
    plot(Ax,By,'b*')
    hold off
end
%% If parallelogram
if any(con3)
    a = [x21 x31; y21 y31];
    b = [Px - x1; Py - y1];

    st = b\a;

    sp = st(:,1);
    tp = st(:,2);

    s = sp(con3);
    t = tp(con3);
end
V = P1.*(1-s).*(1-t) + P2.*s.*(1-t) + P3.*(1-s).*t + P.*s.*t;