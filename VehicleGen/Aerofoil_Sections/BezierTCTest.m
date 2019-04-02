% function [foilData,cLine,tLine] = BezierTC(controlPointsArray,nParts,n,nPop,disc,uDisc)
%Bezier Curves
%Script to plot an aerofoil section using Bezier Curves
%initialise symbolic variable

% if nargin >= 6
%     
%     lDisc = disc;
%     
% elseif nargin >= 5
%     
%     if length(disc) > 1
%         
%         lDisc = disc;
%         
%     elseif disc < 1
%         
%         lDisc = (0:disc:1)';
%     else
%         lDisc = (0:1/(disc - 1):1)';
%     end
%     
%     uDisc = flipud(lDisc);
% else

xPanels = 100;
X = (0:xPanels)';

distribution = "Cosine";

switch distribution
    
    case "Linear"
        
        xNorm = X/max(X);
        
    case "Cosine"
        
        xNorm = 0.5*(1-cos((X*pi)/max(X)));
        
    case "HalfCosine"
        
        xNorm = 1-cos((X*(pi/2))/max(X));
end

lDisc = xNorm;
uDisc = flipud(lDisc);

nPop = 1000;
nParts = 1;

varMin = [1, 0.7,   0.5,    0.3,    0.1,    0,      0;  % xc
    1,      0.7,    0.5,    0.3,    0.1,    0,      0;  % xt
    -0.03, -0.1,   -0.1,   -0.05,  -0.05,   0,      0;  % zc
    0,     -0.05,  -0.05,  -0.05,  -0.05,   0.01,   0]'; % zt

varMax = [1, 0.9,   0.7,    0.5,    0.3,    0,      0;  % xc
    1,      0.9,    0.7,    0.5,    0.3,    0,      0;  % xt
    0.03,   0.1,    0.2,    0.2,    0.2,    0.01,   0;  % zc
    0.03,   0.1,    0.2,    0.2,    0.1,    0.03,   0]'; % zt

n = 6;
nVar = numel(varMin);

varMinMat = repmat(varMin(:)',nPop,1);
varMaxMat = repmat(varMax(:)',nPop,1);

varSize = [nPop nVar];

controlPointsArray = unifrnd(varMinMat,varMaxMat,varSize);

[~,total] = size(controlPointsArray);
array = 1:total;

nSecs = 1;

totPerSec = total/nSecs;

mat = reshape(array,totPerSec,nSecs);

[foilData,cLine,tLine] = deal(cell(nPop,nSecs));

nBez = n + 1;

k = 0:n;
f = factorial(n)./(factorial(k).*factorial(n-k));
s = (0:0.01:1)';

% Only constrain angles past quarter chrod
q3Chord = ceil(length(s)/4);

bezEq = f.*((1-s).^(n-k)).*(s.^k);

for i = 1:nSecs
    
    split = reshape(mat(:,i),nBez,[]);
    
    for j = 1:nPop
        
        % camber x coordinates of control points
        xc = controlPointsArray(j,split(:,1));
        % thickness x coordinates of control points
        xt = controlPointsArray(j,split(:,2));
        % camber z coordinates of control points
        zc = controlPointsArray(j,split(:,3));
        % thickness z coordinates of control points
        zt = controlPointsArray(j,split(:,4));
        %check number of x co-ordinates is equal to number of y co-ordinates
        
        bc = bezEq.*xc;
        cc = bezEq.*zc;
        bt = bezEq.*xt;
        ct = bezEq.*zt;
        
        Pxc = flipud(sum(bc,2));
        Pzc = flipud(sum(cc,2));
        Pxt = flipud(sum(bt,2));
        Pzt = flipud(sum(ct,2));
        
        %% Camber line feasibility test
        xdiff = diff(Pxc);
        zdiff = diff(Pzc);
        
        xdiff(end+1) = xdiff(end);
        zdiff(end+1) = zdiff(end);
        
%         xdiff = [-diff(Pxc); Pxc(end-1) + Pxc(end)];
%         zdiff = [-diff(Pzc); Pzc(end-1) - Pzc(end)];

        theta = atan(abs(zdiff)./xdiff);
        
        theta(isnan(theta)) = 0;
        
        con = zdiff < 0;
        theta(con) = -theta(con);
        
        if any(abs(theta(q3Chord:end)) > 60 * pi/180)
            
            flag(i) = true;
            continue
        end
        
%         theta(end) = 0;
        
        %% Thickness line feasibility test
        xdiff = diff(Pxt);
        zdiff = diff(Pzt);
        
        phi = atan(abs(zdiff)./xdiff);
        
        phi(isnan(phi)) = 0;
        
        con = zdiff < 0;
        phi(con) = -phi(con);
        
        if any(abs(phi(q3Chord:end)) > 50 * pi/180) || any(Pzt(2:end-1) <= 0)
            
            flag(i) = true;
            continue
        end
        
        xu = Pxc - Pzt .* sin(theta);% - t .* sin(theta);
        zu = Pzc + Pzt .* cos(theta);% + t .* sin(theta);
        
        xl = Pxc + Pzt .* sin(theta);% - t .* sin(theta);
        zl = Pzc - Pzt .* cos(theta);
        
        con = [false; xu(2:end) == 0];
        xu(con) = [];
        zu(con) = [];
        
        con = [false; xl(2:end) == 0];
        xl(con) = [];
        zl(con) = [];

        con = [xu(1:end-1) == 1; false];
        xu(con) = [];
        zu(con) = [];
        
        con = [xl(1:end-1) == 1; false];
        xl(con) = [];
        zl(con) = [];
        
        Pzu = interp1(xu, zu, uDisc, 'pchip');
        Pzl = interp1(xl, zl, lDisc, 'pchip');
        
        trailingThick = abs(zu(1) - zl(1));
        
        if trailingThick >= 0.03
        
            continue
        end
        
        [maxThickness,~] = max(abs(Pzu - Pzl));
        
        if maxThickness < 0.1
            
            continue
        end
        
%         chordMaxThickness = uDisc(ID);
        
        % plot aerofoil
        figure
        hold on
        axis equal
        plot(Pxc,Pzc,'r');
        plot(Pxt,Pzt,'b');
        plot(uDisc,Pzu,'k');
        plot(lDisc,Pzl,'k');
        plot(xc,zc,'rx');
        plot(xt,zt,'bx');
        xlabel('x/c');
        ylabel('y/c');
        title('Aerofoil Section');
        legend({'Camber','Thickness','Upper','Lower','Camber Control Points','Thickness Control Points'},'Location','northeast');
        hold off
        
        % 1:end-1 avoids two leading edge points
        foilData{j,i} = [uDisc(1:end-1), Pzu(1:end-1); lDisc, Pzl];
        
        cLine{j,i} = [Pxc Pzc];
        tLine{j,i} = [Pxt Pzt];
    end
end