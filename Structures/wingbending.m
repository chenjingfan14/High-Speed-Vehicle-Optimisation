function M = wingbending(Cp,part,flow)

xyAngle = flow.planeAngles(1);
xzAngle = flow.planeAngles(2);
yzAngle = flow.planeAngles(3);

Uinf = flow.Uinf;
rho = flow.rho;

area = part.area;
centre = part.centre;
unitNorm = part.unitNorm;

[dim1,dim2] = size(centre);

dim2 = dim2/3;

ClAref = -((Cp .* area .* unitNorm(:,:,1)) * sin(xyAngle)) -...
    ((Cp .* area .* unitNorm(:,:,2)) * sin(xzAngle)) -...
    ((Cp .* area .* unitNorm(:,:,3)) * sin(yzAngle));

L = 0.5*rho.*ClAref.*Uinf.^2;

% Offset to normalise z component
centre(:,:,3) = centre(:,:,3) - centre(1,1,3);

yzMag = (centre(:,:,2).^2 + centre(:,:,3).^2).^0.5;

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
    