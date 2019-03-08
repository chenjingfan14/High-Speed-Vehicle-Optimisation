function part = cornervelocities(part,flow)

% clear all
% close all
%
% flow = flowparameters();
% points = plategen();

Vinf = flow.U;

xyz = part.Points;

centre = part.centre;

sym(:,:,1) = centre(:,[1 end],1);
sym(:,:,2) = -centre(:,[1 end],2);
sym(:,:,3) = centre(:,[1 end],3);

centre = [sym(:,1,:), centre, sym(:,end,:)];

area = part.area;

% Adding repeated entries at all edges
pa = xyz([1 1:end-1],:,:) - xyz;
pb = xyz(:,[1 1:end-1],:) - xyz;
pc = xyz([2:end end],:,:) - xyz;
pd = xyz(:,[2:end end],:) - xyz;

pa = centre([1:end end],2:end,:) - xyz;
pb = centre([1 1:end],2:end,:) - xyz;
pc = centre([1 1:end],1:end-1,:) - xyz;
pd = centre([1:end end],1:end-1,:) - xyz;

norm1 = crossmat(pa, pb);
norm2 = crossmat(pb, pc);
norm3 = crossmat(pc, pd);
norm4 = crossmat(pd, pa);

magNorm1 = magmat(norm1);
magNorm2 = magmat(norm2);
magNorm3 = magmat(norm3);
magNorm4 = magmat(norm4);

unitNorm1 = norm1./magNorm1;
unitNorm2 = norm2./magNorm2;
unitNorm3 = norm3./magNorm3;
unitNorm4 = norm4./magNorm4;

con1 = isnan(unitNorm1) & ~isnan(xyz(:,:,1));
con2 = isnan(unitNorm2) & ~isnan(xyz(:,:,1));
con3 = isnan(unitNorm3) & ~isnan(xyz(:,:,1));
con4 = isnan(unitNorm4) & ~isnan(xyz(:,:,1));

unitNorm1(con1) = 0;
unitNorm2(con2) = 0;
unitNorm3(con3) = 0;
unitNorm4(con4) = 0;

con1 = magNorm1 > 0;
con2 = magNorm2 > 0;
con3 = magNorm3 > 0;
con4 = magNorm4 > 0;

numNorms = con1 + con2 + con3 + con4;
norm = (unitNorm1 + unitNorm2 + unitNorm3 + unitNorm4)./numNorms;
magNorm = magmat(norm);

unitNorm = norm./magNorm;

%%
% Nodes that do not have enough adjacent nodes to provide normals
real = ~isnan(xyz(:,:,1));

con = (isnan(magNorm) | magNorm == 0) & real;

[row,col,~] = size(xyz);

array = 1:row;

for i = col:-1:1
    
    fix = array(con(:,i));
    
    first = find(real(:,i),1,'first');
    firstPanel(i) = first;
    
    if isempty(fix)
        
        continue
    end
    
    last = find(real(:,i),1,'last');
    
    % Symmetry plane: Create nodes on negative y-axis
    if any(fix == first)
        
        unitNorm(first,i,1) = -1;
        unitNorm(first,i,2) = 0;
        unitNorm(first,i,3) = 0;
    end
    
    if any(fix == last)
        
        unitNorm(last,i,1) = 1;
        unitNorm(last,i,2) = 0;
        unitNorm(last,i,3) = 0;
    end
end

%% Plotting proof
nx = unitNorm(:,:,1);
ny = unitNorm(:,:,2);
nz = unitNorm(:,:,3);

x = xyz(:,:,1);
y = xyz(:,:,2);
z = xyz(:,:,3);

plotter(xyz)

figure(1)
hold on
axis equal
for i=1:numel(x)
    
    plot3([x(i) x(i) + nx(i)],[y(i) y(i) + ny(i)],[z(i) z(i) + nz(i)],'r');
end
hold off

%% Centre point velocities

T = crossmat(unitNorm, permute(Vinf,[3 1 2]));
Vc = crossmat(T, unitNorm);

part.FirstPanel = firstPanel;
part.Vc = Vc;