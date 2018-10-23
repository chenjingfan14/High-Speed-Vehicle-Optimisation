function s=Overpressure(s,Nose)

global Minf gamma

d=Nose.Radius*2;
Cd=0.91; A=0.067; B=0.44; % J=1;

[Cp1,row]=max([s(1).Cp;s(2).Cp]);
[spCp,col]=max(Cp1); row=row(col);
if row<size(s(1).Cp,1)
    part=1;
else
    part=2; row=row-size(s(1).Cp,1);
end
cols=col:col+1;
con=s(part).y(row,cols)==0;
if isequal(con,[1,1])
    [xsp,i]=min(s(part).x(row,cols(con)));
    i=cols(i);
else
    i=cols(con);
    xsp=s(part).x(row,i);
end
ysp=s(part).y(row,i); zsp=s(part).z(row,i);

for ii=1:length(s)
    part=s(ii);
    cent=part.centre;
    area=part.area;
    [dim1,dim2]=size(cent);
    for i=1:dim1
        for a=1:dim2/3
            j=onetothree(a);
            p=cent(i,j);
            mag=sqrt((xsp-p(1))^2+(ysp-p(2))^2+(zsp-p(3))^2);
            pratio=(A*(Minf^2)*((Cd^0.5)/(mag/d)))+B;
            s(ii).Cpinc(i,a)=(2/(gamma*Minf^2))*(pratio-1);
        end
    end
    s(ii).Cp=s(ii).Cp+s(ii).Cpinc;
end