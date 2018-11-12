function pressureplot(varargin)

figure
hold on
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

for ii=1:length(varargin)
    part=varargin{ii};
    for jj=1:size(part,1)*size(part,2)
        Cp(:,jj)=[max(part(jj).Cp(:)),min(part(jj).Cp(:))]';
    end
    Cpm(:,ii)=[max(Cp(1,:)),min(Cp(2,:))]';
    clear Cp
end
Cpmax=max(Cpm(1,:));
Cpmin=min(Cpm(2,:));
range=Cpmax-Cpmin;
r=[1,0,0];
g=[0,1,0];
b=[0,0,1];

for ii=1:length(varargin)
    part=varargin{ii};
    for jj=1:size(part,1)*size(part,2)
        a=part(jj).centre;
        Cp=part(jj).Cp;
        xyz=part(jj).xyz;
        [x,y]=size(a);
        for i=1:x
            for j=1:3:y
                p=[xyz(i:i+1,j:j+2);xyz(i+1:-1:i,j+3:j+5)];
                psym=p; psym(:,2:3:end)=-psym(:,2:3:end);
                
                if Cp(i,(j+2)/3) <= range/4
                    c1=0;
                    c2=interp1([Cpmin,range/4],[0,1],Cp(i,(j+2)/3));
                    c3=1;
                elseif Cp(i,(j+2)/3) <= range/2
                    c1=0;
                    c2=1;
                    c3=interp1([range/4,range/2],[1,0],Cp(i,(j+2)/3));
                elseif Cp(i,(j+2)/3) <= range*3/4
                    c1=interp1([range/2,range*3/4],[0,1],Cp(i,(j+2)/3));
                    c2=1;
                    c3=0;
                else
                    c1=1;
                    c2=interp1([range*3/4,Cpmax],[1,0],Cp(i,(j+2)/3));
                    c3=0;
                end
                
                colour=[c1,c2,c3];
                fill3(p(:,1),p(:,2),p(:,3),colour,'EdgeColor','none');
                fill3(psym(:,1),psym(:,2),psym(:,3),colour,'EdgeColor','none');
                
            end
        end
    end
end
hold off
end