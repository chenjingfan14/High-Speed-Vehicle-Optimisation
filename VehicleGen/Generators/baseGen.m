function bodybase=baseGen(body,r)

bodybase.Conical = true;

[dim1,dim2] = size(r);

x = repmat(body.x(end),size(r,2),1);

point = zeros(dim1,dim2);

bodybase.Points.Name = "base";
bodybase.Points.xyz = [x,point',x,r'];

bodybase.Points = pointstoxyz(bodybase.Points);