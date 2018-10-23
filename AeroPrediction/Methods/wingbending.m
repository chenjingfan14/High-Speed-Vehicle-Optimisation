function M = wingbending(Cp,points,flow)

xyAngle = flow.planeAngles(1);
xzAngle = flow.planeAngles(2);
yzAngle = flow.planeAngles(3);

Uinf = flow.Uinf;
rho = flow.rho;

area = points.area;
centre = points.centre;
Norm = points.unitNorm;

[dim1,dim2] = size(centre);

X = 1:3:dim2;
Y = X + 1;
Z = Y + 1;

dim2 = dim2/3;

y = centre(:,Y);
z = centre(:,Z);

unitNx = Norm(:,X);
unitNy = Norm(:,Y);
unitNz = Norm(:,Z);

ClAref = -((Cp.*area.*unitNx)*sin(xyAngle)) - ((Cp.*area.*unitNy)*sin(xzAngle)) - ((Cp.*area.*unitNz)*sin(yzAngle));
L = 0.5*rho.*ClAref.*Uinf.^2;

% Offset to normalise z component
z = z - z(1);

yzMag = (y.^2 + z.^2).^0.5;

yzChordAvg = sum(yzMag,1)/dim1;
FChordAvg = sum(L,1);


array = 1:dim2;
mid = array(ceil(end/2));
up = 1:mid-1;
lo = dim2:-1:mid+1;

yzUpper = yzChordAvg(up);
yzLower = yzChordAvg(lo);

yz = [(yzUpper + yzLower)/2, yzChordAvg(mid)];

Fupper = FChordAvg(up);
Flower = FChordAvg(lo);

F = [sum([Fupper;Flower],1), FChordAvg(mid)];

dummy = zeros(mid,1);
S = dummy;
M = dummy;

for i = mid-1:-1:1
    S(i) = S(i+1) - (yz(i+1) - yz(i)).*(F(i+1) + F(i))/2;
    M(i) = M(i+1) - (yz(i+1) - yz(i)).*(S(i+1) + S(i))/2;
end

% yznorm = yz/max(yz);

% figure
% hold on
% plot(yznorm,F,'r')
% plot(yznorm,M,'b')
% hold off

M = sum(M);
    