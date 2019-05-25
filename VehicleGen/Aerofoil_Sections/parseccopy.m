%% PARSEC
%{
This function determines a=[a1, a2, ...an] to solve the airfoil polynomial.
Zn=an(p)*X^(n-1/2), where n is the number of coordinates for the upper or
lower surface.
Input is a vector of PARSEC parameters p=[p1, p2, ...pn] where
p1=rle
p2=Xup
p3=Yup
p4=YXXup
p5=Xlow
p6=Ylow
p7=YXXlow
p8=yte
p9=delta yte (t.e. thickness)
p10=alpha te
p11=beta te
%}

function [zup,zlo] = parsec(parPos,varArray,n,nPop)

x = 0:0.01:1;

rleu = parPos(:,varArray == "rleu");
% rlel = parPos(:,varArray == "rlel");
xup = parPos(:,varArray == "xup");
zup = parPos(:,varArray == "zup");
zxxup = parPos(:,varArray == "zxxup");
xlo = parPos(:,varArray == "xlo");
zlo = parPos(:,varArray == "zlo");
zxxlo = parPos(:,varArray == "zxxlo");
zte = parPos(:,varArray == "zte");
dzte = parPos(:,varArray == "dzte");
ate = parPos(:,varArray == "ate");
bte = parPos(:,varArray == "bte");

for i = 1:nPop * n
    
    p(1) = rleu(i);
    % rlel = parPos(:,varArray == "rlel");
    p(2) = xup(i);
    p(3) = zup(i);
    p(4) = zxxup(i);
    p(5) = xlo(i);
    p(6) = zlo(i);
    p(7) = zxxlo(i);
    p(8) = zte(i);
    p(9) = dzte(i);
    p(10) = ate(i);
    p(11) = bte(i);
    
    % Define matricies
    c1=[1,1,1,1,1,1];
    c2=[p(2)^(1/2),p(2)^(3/2),p(2)^(5/2),p(2)^(7/2),p(2)^(9/2),p(2)^(11/2)];
    c3=[1/2, 3/2, 5/2, 7/2, 9/2, 11/2];
    c4=[(1/2)*p(2)^(-1/2), (3/2)*p(2)^(1/2),(5/2)*p(2)^(3/2),(7/2)...
        *p(2)^(5/2),(9/2)*p(2)^(7/2),(11/2)*p(2)^(9/2)];
    c5=[(-1/4)*p(2)^(-3/2),(3/4)*p(2)^(-1/2),(15/4)*p(2)^(1/2),(35/4)...
        *p(2)^(3/2),(53/4)*p(2)^(5/2),(99/4)*p(2)^(7/2)];
    c6=[1,0,0,0,0,0];
    Cup=[c1; c2; c3; c4; c5; c6];
    c7=[1,1,1,1,1,1];
    c8=[p(5)^(1/2),p(5)^(3/2),p(5)^(5/2),p(5)^(7/2),p(5)^(9/2),p(5)^(11/2)];
    c9=[1/2, 3/2, 5/2, 7/2, 9/2, 11/2];
    c10=[(1/2)*p(5)^(-1/2), (3/2)*p(5)^(1/2),(5/2)*p(5)^(3/2),(7/2)...
        *p(5)^(5/2),(9/2)*p(5)^(7/2),(11/2)*p(5)^(9/2)];
    c11=[(-1/4)*p(5)^(-3/2),(3/4)*p(5)^(-1/2),(15/4)*p(5)^(1/2),(35/4)...
        *p(5)^(3/2),(53/4)*p(5)^(5/2),(99/4)*p(5)^(7/2)];
    c12=[1,0,0,0,0,0];
    Clo=[c7; c8; c9; c10; c11; c12];
    bup=[p(8)+p(9)/2;p(3);tand(-p(10)+p(11)/2);0;p(4);(sqrt(2*p(1)))];
    blo=[p(8)-p(9)/2;p(6);tand(p(10)+p(11)/2);0;p(7);-(sqrt(2*p(1)))];
    % Solve system of equations: C x a=b
    aup=linsolve(Cup,bup);
    alo=linsolve(Clo,blo);
    % Format a
    
    aup(:,1)=aup;
    alo(7:12,1)=alo;
    
    a1 = aup(1);
    a2 = aup(2);
    a3 = aup(3);
    a4 = aup(4);
    a5 = aup(5);
    a6 = aup(6);
    
    zup = a6.*x.^5.5 + a5.*x.^4.5 + a4.*x.^3.5 + a3.*x.^2.5 + a2.*x.^1.5 + a1.*x.^0.5;
    
    a1 = alo(1);
    a2 = alo(2);
    a3 = alo(3);
    a4 = alo(4);
    a5 = alo(5);
    a6 = alo(6);
    
    zlo = a6.*x.^5.5 + a5.*x.^4.5 + a4.*x.^3.5 + a3.*x.^2.5 + a2.*x.^1.5 + a1.*x.^0.5;
    
    plot(x,zup,'k',x,zlo,'k')
    
end