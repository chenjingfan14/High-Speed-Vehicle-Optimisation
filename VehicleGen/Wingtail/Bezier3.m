function foilData = Bezier3(controlPointsArray,nParts,n,nPop,disc,uDisc)
%Bezier Curves
%Script to plot an aerofoil section using Bezier Curves
%initialise symbolic variable

if nargin >= 6
    
    lDisc = disc;
    
elseif nargin >= 5
    
    if length(disc) > 1
        
        lDisc = disc;
        
    elseif disc < 1
        
        lDisc = (0:disc:1)';
    else
        lDisc = (0:1/(disc - 1):1)';
    end
    
    uDisc = flipud(lDisc);
else
    lDisc = (0:1/100:1)';
    uDisc = flipud(lDisc);
end

[~,total] = size(controlPointsArray);
array = 1:total;

nSecs = nParts + 1;

totPerSec = total/nSecs;

mat = reshape(array,totPerSec,nSecs);

foilData = cell(nPop,nSecs);

nBez = n + 1;

k = 0:n;
f = factorial(n)./(factorial(k).*factorial(n-k));
s = (0:0.01:1)';

bezEq = f.*((1-s).^(n-k)).*(s.^k);

for i = 1:nSecs
    
    split = reshape(mat(:,i),nBez,[]);
    
    for j = 1:nPop
        %set i for required number of aerofoil configurations
        % upper x coordinates of control points
        
        xu = controlPointsArray(j,split(:,1));
        % lower x coordinates of control points
        xl = controlPointsArray(j,split(:,2));
        % upper y coordinates of control points
        yu = controlPointsArray(j,split(:,3));
        % lower y coordinates of control points
        yl = controlPointsArray(j,split(:,4));
        %check number of x co-ordinates is equal to number of y co-ordinates
        
        bu = bezEq.*xu;
        cu = bezEq.*yu;
        bl = bezEq.*xl;
        cl = bezEq.*yl;
        
        Pxu = sum(bu,2);
        Pzu = sum(cu,2);
        Pxl = sum(bl,2);
        Pzl = sum(cl,2);
        
        Pzu = interp1(Pxu, Pzu, uDisc);
        Pzl = interp1(Pxl, Pzl, lDisc);
        
        % plot aerofoil
%         figure
%         hold on
%         axis equal
%         plot(uDisc,Pzu,'r');
%         plot(lDisc,Pzl,'b');
%         plot(xu,yu,'rx');
%         plot(xl,yl,'bx');
%         xlabel('x/c');
%         ylabel('y/c');
%         title('Aerofoil Section');
%         axis([0 1 -0.1 0.1]);
%         legend({'Upper','Lower','Upper Control Points','Lower Control Points'},'Location','northeast');
%         hold off
        
        % 1:end-1 avoids two leading edge points
        foilData{j,i} = [uDisc(1:end-1), Pzu(1:end-1); lDisc, Pzl];
        
    end
end