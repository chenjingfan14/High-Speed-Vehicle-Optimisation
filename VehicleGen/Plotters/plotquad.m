function plotquad(a,b,c,d)

if nargin == 1

    x = a(:,:,1);
    y = a(:,:,2);
    z = a(:,:,3);
    
    x = x([1 2 4 3 1]);
    y = y([1 2 4 3 1]);
    z = z([1 2 4 3 1]);
    
else
    
    x = [a(1) b(1) c(1) d(1) a(1)];
    y = [a(2) b(2) c(2) d(2) a(2)];
    z = [a(3) b(3) c(3) d(3) a(3)];
end

plot3(x,y,z)