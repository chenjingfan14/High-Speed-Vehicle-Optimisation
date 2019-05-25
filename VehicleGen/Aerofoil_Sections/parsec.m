function [aup,alo] = parsec(parPos,varArray,n,nPop)
%{
This function determines a=[a1, a2, ...an] to solve the airfoil polynomial.
Zn=an(p)*X.^(n-1/2), where n is the number of coordinates for the upper or
lower surface. 
Input is a vector of PARSEC parameters p=[p1, p2, ...pn] where
p1 = rleup
p2 = rlelo
p3 = xup
p4 = zup
p5 = zxxup
p6 = xlo
p7 = zlo
p8 = zxxlo
p9 = zte
p10 = delta zte (t.e. thickness)
p11 = alpha te
p12 = beta te
%}

rleu = permute(parPos(:,varArray == "rleu"), [3 1 2]);
rlel = permute(parPos(:,varArray == "rlel"), [3 1 2]);
xup = permute(parPos(:,varArray == "xup"), [1 3 2]);
zup = permute(parPos(:,varArray == "zup"), [3 1 2]);
zxxup = permute(parPos(:,varArray == "zxxup"), [3 1 2]);
xlo = permute(parPos(:,varArray == "xlo"), [1 3 2]);
zlo = permute(parPos(:,varArray == "zlo"), [3 1 2]);
zxxlo = permute(parPos(:,varArray == "zxxlo"), [3 1 2]);
zte = permute(parPos(:,varArray == "zte"), [3 1 2]);
dzte = permute(parPos(:,varArray == "dzte"), [3 1 2]);
ate = permute(parPos(:,varArray == "ate"), [3 1 2]);
bte = permute(parPos(:,varArray == "bte"), [3 1 2]);

% Define matricies

c2 = [xup.^(1/2), xup.^(3/2), xup.^(5/2), xup.^(7/2), xup.^(9/2), xup.^(11/2)];

c3 = [1/2, 3/2, 5/2, 7/2, 9/2, 11/2];

c4 = [(1/2)*xup.^(-1/2), (3/2)*xup.^(1/2), (5/2)*xup.^(3/2),(7/2)...
    *xup.^(5/2), (9/2)*xup.^(7/2), (11/2)*xup.^(9/2)];

c5 = [(-1/4)*xup.^(-3/2), (3/4)*xup.^(-1/2), (15/4)*xup.^(1/2), (35/4)...
    *xup.^(3/2), (53/4)*xup.^(5/2), (99/4)*xup.^(7/2)];

c1 = [1,1,1,1,1,1];

c6 = [1,0,0,0,0,0];

c7 = [1,1,1,1,1,1];

c8 = [xlo.^(1/2), xlo.^(3/2), xlo.^(5/2), xlo.^(7/2), xlo.^(9/2), xlo.^(11/2)];

c9 = [1/2, 3/2, 5/2, 7/2, 9/2, 11/2];

c10 = [(1/2)*xlo.^(-1/2), (3/2)*xlo.^(1/2), (5/2)*xlo.^(3/2), (7/2)...
    *xlo.^(5/2), (9/2)*xlo.^(7/2), (11/2)*xlo.^(9/2)];

c11 = [(-1/4)*xlo.^(-3/2), (3/4)*xlo.^(-1/2), (15/4)*xlo.^(1/2), (35/4)...
    *xlo.^(3/2), (53/4)*xlo.^(5/2), (99/4)*xlo.^(7/2)];

c12 = [1,0,0,0,0,0];

bup = [zte + dzte; zup; tand(ate + bte/2); zeros(size(zup)); zxxup; (2*rleu).^0.5];
blo = [zte; zlo; tand(ate - bte/2); zeros(size(zup)); zxxlo; (2*rlel).^0.5];

for i = nPop:-1:1
    
    for j = n + 1:-1:1
    
        cup = [c1; c2(i,:,j); c3; c4(i,:,j); c5(i,:,j); c6];

        clo = [c7; c8(i,:,j); c9; c10(i,:,j); c11(i,:,j); c12];

        % Solve system of equations: C x a=b
        aup(i,:,j) = linsolve(cup,bup(:,i,j));
        alo(i,:,j) = linsolve(clo,blo(:,i,j));
    end
end

a1 = aup(:,1,:);
a2 = aup(:,2,:);
a3 = aup(:,3,:);
a4 = aup(:,4,:);
a5 = aup(:,5,:);
a6 = aup(:,6,:);

x = 0:0.01:1;

z = a6.*x.^5.5 + a5.*x.^4.5 + a4.*x.^3.5 + a3.*x.^2.5 + a2.*x.^1.5 + a1.*x.^0.5;
