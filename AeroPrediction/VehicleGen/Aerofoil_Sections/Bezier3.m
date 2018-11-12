function foilData = Bezier3(controlPointsArray,nParts,nBez,nPop)
%Bezier Curves
%Script to plot an aerofoil section using Bezier Curves
%initialise symbolic variable

[~,total] = size(controlPointsArray);
array = 1:total;

nSecs = nParts + 1;

totPerSec = total/nSecs;
totPerVar = total/nBez;
BezPerVar = totPerSec/totPerVar;

mat = reshape(array,totPerSec,nSecs);

foilData = cell(nPop,nSecs);

for i = 1:nSecs
    
    split = reshape(mat(:,i),BezPerVar,totPerVar);
    
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
        
        n = BezPerVar - 1;
        
        k = 0:n;
        f = factorial(n)./(factorial(k).*factorial(n-k));
        
        s = (0:0.01:1)';
        
        bu = f.*((1-s).^(n-k)).*(s.^k).*xu;
        cu = f.*((1-s).^(n-k)).*(s.^k).*yu;
        bl = f.*((1-s).^(n-k)).*(s.^k).*xl;
        cl = f.*((1-s).^(n-k)).*(s.^k).*yl;
        
        Pxu = sum(bu,2);
        Pzu = sum(cu,2);
        Pxl = sum(bl,2);
        Pzl = sum(cl,2);
        
        % plot aerofoil
%         figure
%         hold on
%         axis equal
%         plot(Pxu,Pzu,'r');
%         plot(Pxl,Pzl,'b');
%         plot(xu,yu,'rx');
%         plot(xl,yl,'bx');
%         xlabel('x/c');
%         ylabel('y/c');
%         title('Aerofoil Section');
%         axis([0 1 -0.1 0.1]);
%         legend({'Upper','Lower','Upper Control Points','Lower Control Points'},'Location','northeast');
%         hold off
        
        foilData{j,i} = [Pxu, Pzu; Pxl, Pzl];
        
    end
end