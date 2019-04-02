% function foilData = Bezier3(controlPointsArray,nParts,nBez,nPop)

%Bezier Curves
%Script to plot an aerofoil section using Bezier Curves

nPop = 100;
nParts = 1;

varMin = [0.13 0.08 0.28 0.13 0.43 0.58 0.76 0.03 0.01 0.18 0.06 0.38 0.58 0.78 0];
varMax = [0.22 0.17 0.37 0.22 0.52 0.67 0.85 0.12 0.10 0.27 0.15 0.47 0.67 0.87 0.09];

nVar = length(varMin);
nBez = 7;

varMinMat = repmat(varMin,nPop,1);
varMaxMat = repmat(varMax,nPop,1);

varSize = [nPop nVar];

controlPointsArray = unifrnd(varMinMat,varMaxMat,varSize);
controlPointsArray = varMax;

[~,total] = size(controlPointsArray);
array = 1:total;

% nSecs = nParts + 1;
nSecs = 1;

totPerSec = total/nSecs;

mat = reshape(array,totPerSec,nSecs);

foilData = cell(nPop,nSecs);

s = (0:1/99:1)';
n = nBez - 1;
        
k = 0:n;
f = factorial(n)./(factorial(k).*factorial(n-k));

for i = 1:nSecs
    
    split = mat(:,i);
    
    for j = 1:nPop
        %set i for required number of aerofoil configurations
        % upper x coordinates of control points
        
        x = [0 controlPointsArray(j,split([1 3 5 6 7])) 1];
        y = [0 controlPointsArray(j,split([2 4 4 4 8])) 0];
        xt = [0 0 controlPointsArray(j,split([10 12 13 14])) 1];
        t = [0 controlPointsArray(j,split([9 11 11 11 15])) 0];
        %check number of x co-ordinates is equal to number of y co-ordinates
        
        b = f.*((1-s).^(n-k)).*(s.^k).*x;
        c = f.*((1-s).^(n-k)).*(s.^k).*y;
        
        Px = sum(b,2);
        Pz = sum(c,2);
        
        b = f.*((1-s).^(n-k)).*(s.^k).*xt;
        c = f.*((1-s).^(n-k)).*(s.^k).*t;
        
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