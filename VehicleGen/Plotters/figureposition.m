function pos = figureposition()

h = findobj('type','figure');
num = length(h)+1;

num6 = ceil(num/6);

if num6 > 1
    num = num - 6*(num6-1);
end

n = num./(1:6);
isWhole = rem(n,1) == 0;
lowestWhole = max(n(isWhole));

switch lowestWhole
    case 1
        pos = [0, 560, 640, 440];
    case 3 
        pos = [640, 560, 640, 440];
    case 5
        pos = [1280, 560, 640, 440];
    case 2
        pos = [0, 40, 640, 440];
    case 4
        pos = [640, 40, 640, 440];
    case 6
        pos = [1280, 40, 640, 440];
        
%     case 1
%         pos = [0, -150, 640, 440];
%     case 3 
%         pos = [640, -150, 640, 440];
%     case 5
%         pos = [1280, -150, 640, 440];
%     case 2
%         pos = [0, 380, 640, 440];
%     case 4
%         pos = [640, 380, 640, 440];
%     case 6
%         pos = [1280, 380, 640, 440];
        
end