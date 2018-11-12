function logpressureplot(Assembly,alpha,M)

for ii=1:length(Assembly)
    part=Assembly(ii);
    for jj=1:size(part,1)*size(part,2)
        Cp(:,jj)=[max(part(jj).Cp(:)),min(part(jj).Cp(:))]';
    end
    Cpm(:,ii)=[max(Cp(1,:)),min(Cp(2,:))]';
    clear Cp
end

Cpmax=log(max(Cpm(1,:)));
Cpmin=log(0.0001);
range=-Cpmax+Cpmin;

figure
title(['Log(Cp) at Alpha = ' num2str(alpha) ', M = ' num2str(M) ' (Cpmax = ' num2str(max(Cpm(1,:))) ')']); 
hold on
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

% r=[1,0,0];
% g=[0,1,0];
% b=[0,0,1];

count=1;
for ii=1:length(Assembly)
    part=Assembly(ii);
    for jj=1:size(part,1)*size(part,2)
        a=part(jj).centre;
        Cp=log(part(jj).Cp+1*10^-32);
        xyz=part(jj).xyz;
        [x,y]=size(a);
        for i=1:x
            for j=1:3:y
                p=[xyz(i:i+1,j:j+2);xyz(i+1:-1:i,j+3:j+5)];
                psym=p; psym(:,2:3:end)=-psym(:,2:3:end);
                
                val=Cp(i,(j+2)/3);
                if exp(val)<0.0001
                    val=log(0.0001);
                end
                test(i,(j+2)/3)=val;
                
                if val <= range*3/4
                    c1=0;
                    c2=interp1([Cpmin,range*3/4],[0,1],val);
                    c3=1;
                elseif val <= range/2
                    c1=0;
                    c2=1;
                    c3=interp1([range*3/4,range/2],[1,0],val);
                elseif val <= range/4
                    c1=interp1([range/2,range/4],[0,1],val);
                    c2=1;
                    c3=0;
                else
                    c1=1;
                    c2=interp1([range/4,Cpmax],[1,0],val);
                    c3=0;
                end
                
                colour=[c1,c2,c3];
                fill3(p(:,1),p(:,2),p(:,3),colour,'EdgeColor','none');
                fill3(psym(:,1),psym(:,2),psym(:,3),colour,'EdgeColor','none');
                c(count,:)=colour;
                count=count+1;
            end
        end
    end
end
% Cpmin=exp(Cpmin); Cpmax=exp(Cpmax);
% c=Cpmin:(Cpmax-Cpmin)/8:Cpmax;
% caxis([c(1) c(end)])
% colormap(jet)
% colorbar('Ytick',log(c),'Yticklabel',c)
hold off
end