function points = deflectcontrolsufrace(properties,points,delta,deltaPrev)

if nargin > 3
    delta = delta - deltaPrev;
end

ID = properties.Control.ChordID;
controlSurf = properties.Control.Surf;

x = points.x;
y = points.y;
z = points.z;

di = properties.Dihedral;

[row,~] = size(x);
rows = ID:row;

z(rows,controlSurf) = z(rows,controlSurf) + (x(rows,controlSurf) - x(ID,controlSurf)) * sin(delta) * cos(di);
y(rows,controlSurf) = y(rows,controlSurf) + (x(rows,controlSurf) - x(ID,controlSurf)) * sin(delta) * sin(-di);
x(rows,controlSurf) = (x(rows,controlSurf) - x(ID,controlSurf)) * cos(delta) + x(ID,controlSurf);

points.x = x;
points.y = y;
points.z = z;

points = xyztopoints(points);

end