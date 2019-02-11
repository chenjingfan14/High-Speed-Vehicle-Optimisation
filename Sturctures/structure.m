function structure(wing, part, wingP)
% Initial structural program. See Low-fidelity aerostructural optimization
% of aircraft wings with a simplified wingbox model using OpenAeroStruct

if iscell(wing)
    
    wing = wing{:};
end

%% Initialise Skin
skin = wing.Skin;
nSecs = ceil(size(wingP,2)/2);
di = wing.Dihedral;

skinThick = skin.Thickness;

DoF = 6;

E = 70e9;
v = 0.29;
G = E / (2 * (1 + v));

% Loops through spanwise sections
for i = nSecs:-1:1
    
    skinUpper = skin.Upper{i};
    skinLower = skin.Lower{i};
    skinUpperMean = skin.UpperMean{i};
    skinLowerMean = skin.LowerMean{i};
    
    [row,~] = size(skinUpper);
    j = [2, 2:row-1, row-1];
    
    % This or skinUpperMean etc?
    skinUpOut = squeeze(skinUpper(:,1,:));
    skinLoOut = squeeze(skinLower(:,1,:));
    skinUpIn = squeeze(skinUpper(j,2,:));
    skinLoIn = squeeze(skinLower(j,2,:));
    
    %% Calculate skin mid-line properties
    % Find middle points of skins, use inner second and second last twice
    % instead of first and last to account for spar thickness
    
    midUpper = (skinUpOut + skinUpIn)/2;
    midLower = (skinLoOut + skinLoIn)/2;
    
    a = midUpper(1:end-1,:,:);
    b = midUpper(2:end,:,:);
    c = midLower(2:end,:,:);
    d = midLower(1:end-1,:,:);
    
    skinUpperVec(:,:,i) = b - a;
    skinLowerVec(:,:,i) = c - d;
    
    skinUpperL = ((b(:,1) - a(:,1)).^2 + (b(:,2) - a(:,2)).^2 + (b(:,3) - a(:,3)).^2) .^0.5;
    skinLowerL = ((c(:,1) - d(:,1)).^2 + (c(:,2) - d(:,2)).^2 + (c(:,3) - d(:,3)).^2) .^0.5;
    
    skinLength(:,i,1) = skinUpperL;
    skinLength(:,i,2) = skinLowerL;
    
    %%
    figure(1)
    hold on
    axis equal
    plot3(skinUpper(:,:,1),skinUpper(:,:,2),skinUpper(:,:,3),'r')
    % plot3(xMidUpper,yMidUpper,zMidUpper,'r--')
    plot3(skinLower(:,:,1),skinLower(:,:,2),skinLower(:,:,3),'r')
    % plot3(xMidLower,yMidLower,zMidLower,'r--')
    
end

%% Initiliase Spars (Only works for 4 points per spar)
spars = wing.Spar.Points;

nSpars = numel(spars);

sparThick = wing.Spar.Thickness;

% Loops through chordwise spars
for i = nSpars:-1:1
    
    spar = spars{i};
    
    sparFront = spar(1:2,:,:);
    sparBack = spar(3:4,:,:);
    sparMid = [mean([sparFront(1,:,:); sparBack(1,:,:)]);...
        mean([sparFront(2,:,:); sparBack(2,:,:)])];
    
    %% Calculate spar mid-line properties
    
    sparMidDiff = squeeze(diff(sparMid));
    
    sparMidVec(i,:,:) = sparMidDiff';
end

%% Create third reference nodes


%% Total structural configuration properties

beamVectors = [sparMidVec(1,:,:); skinUpperVec;
    sparMidVec(end,:,:); flipud(skinLowerVec)];

L = squeeze(beamVectors(:,1,:).^2 + beamVectors(:,2,:).^2 + beamVectors(:,3,:).^2).^0.5;

r = [sparThick(1,:); repmat(skinThick, size(skinUpperVec,1), 1);
    sparThick(end,:); repmat(skinThick, size(skinLowerVec,1), 1)];

A = pi * r.^2;
Iz = pi * r.^4 / 4;
Iy = Iz;

J = Iy + Iz;

nBeams = size(L,1);

%% Replace by creating box points earlier?
box = wing.Box;

first = find(box,1,'first');
last = find(box,1,'last');

[boxDim,~] = size(box);

%% CHECK: Points correspond to panel coordinates ahead of panel centre
wingPoints = part.Points;
wingCentre = part.centre;
wingS = part.area;

boxPoints = wingPoints(box,:,:);
boxCentre = wingCentre(box,:,:);
boxP = wingP(box,:,:);
boxS = wingS(box,:);

node1 = boxPoints(:,1:end-1,:);
node2 = boxPoints(:,2:end,:);

node1 = [repmat(node1(1,:,:), first, 1); node1(2:end-1,:,:); repmat(node1(end,:,:), boxDim - last + 1, 1)];
node2 = [repmat(node2(1,:,:), first, 1); node2(2:end-1,:,:); repmat(node2(end,:,:), boxDim - last + 1, 1)];

force = wingP .* wingS;

vec1 = [wingCentre; node1(end,:,:)] - node1;
vec2 = [wingCentre; node2(end,:,:)] - node2;
vec3 = [node1(1,:,:); wingCentre] - node1;
vec4 = [node2(1,:,:); wingCentre] - node2;

% Adding imaginary nodes to depict zero force felt inward of wing
% forceLeft = [zeros(row,1,3), force];
% forceRight = [force, zeros(row,1,3)];
% vec1 = [zeros(row,1,3), vec1];
% vec2 = [vec1, zeros(row,1,3)];

% Begin/end forces repeated for LE/TE forces
% Symmetry plane used to make inboard wing have identical forces from left
% & right
forceLeft1 = force([1:end end],[1 1:end],:);
forceRight1 = force([1:end end],[1:end end],:);

forceLeft2 = force([1 1:end],[1 1:end],:);
forceRight2 = force([1 1:end],[1:end end],:);

% Must invert y-direction for this
symForce1 = vec1(:,1,:);
symForce2 = vec2(:,end,:);

symForce1(:,:,2) = -symForce1(:,:,2);
symForce2(:,:,2) = -symForce2(:,:,2);

vec1sym = [symForce1, vec1];
vec2sym = [vec2, symForce2];
vec3sym = [symForce1, vec3];
vec4sym = [vec4, symForce2];

F = 0.25 * forceLeft1 + 0.25 * forceRight1 +...
    0.25 * forceLeft2 + 0.25 * forceRight2;

crossLeft1(:,:,1) = vec1sym(:,:,2).*forceLeft1(:,:,3) - vec1sym(:,:,3).*forceLeft1(:,:,2);
crossLeft1(:,:,2) = vec1sym(:,:,3).*forceLeft1(:,:,1) - vec1sym(:,:,1).*forceLeft1(:,:,3);
crossLeft1(:,:,3) = vec1sym(:,:,1).*forceLeft1(:,:,2) - vec1sym(:,:,2).*forceLeft1(:,:,1);

crossRight1(:,:,1) = vec2sym(:,:,2).*forceRight1(:,:,3) - vec2sym(:,:,3).*forceRight1(:,:,2);
crossRight1(:,:,2) = vec2sym(:,:,3).*forceRight1(:,:,1) - vec2sym(:,:,1).*forceRight1(:,:,3);
crossRight1(:,:,3) = vec2sym(:,:,1).*forceRight1(:,:,2) - vec2sym(:,:,2).*forceRight1(:,:,1);

crossLeft2(:,:,1) = vec3sym(:,:,2).*forceLeft2(:,:,3) - vec3sym(:,:,3).*forceLeft2(:,:,2);
crossLeft2(:,:,2) = vec3sym(:,:,3).*forceLeft2(:,:,1) - vec3sym(:,:,1).*forceLeft2(:,:,3);
crossLeft2(:,:,3) = vec3sym(:,:,1).*forceLeft2(:,:,2) - vec3sym(:,:,2).*forceLeft2(:,:,1);

crossRight2(:,:,1) = vec4sym(:,:,2).*forceRight2(:,:,3) - vec4sym(:,:,3).*forceRight2(:,:,2);
crossRight2(:,:,2) = vec4sym(:,:,3).*forceRight2(:,:,1) - vec4sym(:,:,1).*forceRight2(:,:,3);
crossRight2(:,:,3) = vec4sym(:,:,1).*forceRight2(:,:,2) - vec4sym(:,:,2).*forceRight2(:,:,1);

M = 0.25 * crossLeft1 + 0.25 * crossRight1 +...
    0.25 * crossLeft2 + 0.25 * crossRight2 ;

%% Sum forces and moments at leading and trailing edge

F = [sum(F(1:first,:,:),1); F(first+1:last-1,:,:); sum(F(last:end,:,:),1)];
M = [sum(M(1:first,:,:),1); M(first+1:last-1,:,:); sum(M(last:end,:,:),1)];

% Rearrange wing spanwise wrap to wingbox chordwisel wrap 
[row,col,~] = size(F);
dim = row*col*3;

halfCol = ceil(col/2);

origID = zeros(row, col, 3);
origID(:) = 1:dim;

transID = [origID(:, 1:halfCol, :); rot90(origID(:, halfCol + 1:end, :),2)];

F = F(transID);
M = M(transID);

% Reshape points, put first node at end to complete circle
boxPoints = [boxPoints(transID); boxPoints(transID(1,:,:))];
node1 = boxPoints(1:end-1,:,:);
node2 = boxPoints(2:end,:,:);

node3 = (node1 + node2)/2;
r1 = node2 - node1;

% No dihedral applied to inboard section, only needs to be altered in y dim
node3(:,1,2) = node3(:,1,2) + 0.5 * L(:,1);
% node3(:,2:end,[2 3]) = node3(:,2:end,[2 3]) + 0.5 * permute(L(:,2:end),[1 3 2]) .* permute(repmat([cos(di); sin(di)], 1, nSecs-1), [3 2 1]);
node3(:,2:end,[2 3]) = node3(:,2:end,[2 3]) + 0.5 * L(:,2:end) .* permute(repmat([cos(di); sin(di)], 1, nSecs-1), [3 2 1]);

r13 = node3 - node1;

r3(:,:,1) = r1(:,:,2).*r13(:,:,3) - r1(:,:,3).*r13(:,:,2);
r3(:,:,2) = r1(:,:,3).*r13(:,:,1) - r1(:,:,1).*r13(:,:,3);
r3(:,:,3) = r1(:,:,1).*r13(:,:,2) - r1(:,:,2).*r13(:,:,1);

r2(:,:,1) = r3(:,:,2).*r1(:,:,3) - r3(:,:,3).*r1(:,:,2);
r2(:,:,2) = r3(:,:,3).*r1(:,:,1) - r3(:,:,1).*r1(:,:,3);
r2(:,:,3) = r3(:,:,1).*r1(:,:,2) - r3(:,:,2).*r1(:,:,1);

r1Mag = (r1(:,:,1).^2 + r1(:,:,2).^2 + r1(:,:,3).^2).^0.5;
r2Mag = (r2(:,:,1).^2 + r2(:,:,2).^2 + r2(:,:,3).^2).^0.5;
r3Mag = (r3(:,:,1).^2 + r3(:,:,2).^2 + r3(:,:,3).^2).^0.5;

lmn1 = r1 ./ r1Mag;
lmn2 = r2 ./ r2Mag;
lmn3 = r3 ./ r3Mag;

AE_L = A .* E ./ L;
E12_L3 = 12 * E ./ L.^3;
E6_L2 = 6 * E ./ L.^2;
GJ_L = G * J ./ L;
E4_L = 4 * E ./ L;
E2_L = 2 * E ./ L;

%% Assemble load vector
xForce = 1:DoF:dim/nSecs*2;
yForce = xForce + 1;
zForce = yForce + 1;
xMoment = zForce + 1;
yMoment = xMoment + 1;
zMoment = yMoment + 1;

[loads,qf] = deal(zeros(dim/nSecs*2,nSecs));

for i = nSecs:-1:1
    
    loads(xForce,i) = F(:,i,1);
    loads(yForce,i) = F(:,i,2);
    loads(zForce,i) = F(:,i,3);
    loads(xMoment,i) = M(:,i,1);
    loads(yMoment,i) = M(:,i,2);
    loads(zMoment,i) = M(:,i,3);
    
%     loads(:) = 1;

    AE_Li = permute(AE_L(:,i), [3 2 1]);
    E12_L3i = permute(E12_L3(:,i), [3 2 1]);
    E6_L2i = permute(E6_L2(:,i), [3 2 1]);
    GJ_Li = permute(GJ_L(:,i), [3 2 1]);
    E4_Li = permute(E4_L(:,i), [3 2 1]);
    E2_Li = permute(E2_L(:,i), [3 2 1]);

    Iyi = permute(Iy(:,i), [3 2 1]);
    Izi = permute(Iz(:,i), [3 2 1]);
    
    lambdaMat = permute([lmn1(:,i,:), lmn2(:,i,:), lmn3(:,i,:)],[2 3 1]);
    
    l1 = lambdaMat(1,1,:);
    l2 = lambdaMat(2,1,:);
    l3 = lambdaMat(3,1,:);
    m1 = lambdaMat(1,2,:);
    m2 = lambdaMat(2,2,:);
    m3 = lambdaMat(3,2,:);
    n1 = lambdaMat(1,3,:);
    n2 = lambdaMat(2,3,:);
    n3 = lambdaMat(3,3,:);
    
    l12 = l1.^2;
    l22 = l2.^2;
    l32 = l3.^2;
    m12 = m1.^2;
    m22 = m2.^2;
    m32 = m3.^2;
    n12 = n1.^2;
    n22 = n2.^2;
    n32 = n3.^2;
    
    k1 = [(AE_Li .* l12) + E12_L3i.*(Izi .* l22 + Iyi .* l32), (AE_Li .* l1 .* m1) + E12_L3i.*(Izi .* l2 .* m2 + Iyi .* l3 .* m3), (AE_Li .* l1 .* n1) + E12_L3i.*(Izi .* l2 .* n2 + Iyi .* l3 .* n3);...
        (AE_Li .* l1 .* m1) + E12_L3i.*(Izi .* l2 .* m2 + Iyi .* l3 .* m3), (AE_Li .* m12) + E12_L3i.*(Izi .* m22 + Iyi .* m32), (AE_Li .* m1 .* n1) + E12_L3i.*(Izi .* m2 .* n2 + Iyi .* m3 .* n3);...
        (AE_Li .* l1 .* n1) + E12_L3i.*(Izi .* l2 .* n2 + Iyi .* l3 .* n3), (AE_Li .* m1 .* n1) + E12_L3i.*(Izi .* m2 .* n2 + Iyi .* m3 .* n3), (AE_Li .* n12) + E12_L3i.*(Izi .* n22 + Iyi .* n32)];
    
    k2 = [E6_L2i.*(Izi - Iyi) .* l2 .* l3, E6_L2i.*(Izi .* l2 .* m3 - Iyi .* l3 .* m2), E6_L2i.*(Izi .* l2 .* n3 - Iyi .* l3 .* n2);...
        E6_L2i.*(Izi .* l3 .* m2 - Iyi .* l2 .* m3), E6_L2i.*(Izi - Iyi) .* m2 .* m3, E6_L2i.*(Izi .* m2 .* n3 - Iyi .* m3 .* n2);...
        E6_L2i.*(Izi .* l3 .* n2 - Iyi .* l2 .* n3), E6_L2i.*(Izi .* m3 .* n2 - Iyi .* m2 .* n3), E6_L2i.*(Izi - Iyi) .* n2 .* n3];
    
    k3 = [(GJ_Li .* l12) + E4_Li.*(Iyi .* l22 + Izi .* l32), (GJ_Li .* l1 .* m1) + E4_Li.*(Iyi .* l2 .* m2 + Izi .* l3 .* m3), (GJ_Li .* l1 .* n1) + E4_Li.*(Iyi .* l2 .* n2 + Izi .* l3 .* n3);...
        (GJ_Li .* l1 .* m1) + E4_Li.*(Iyi .* l2 .* m2 + Izi .* l3 .* m3), (GJ_Li .* m12) + E4_Li.*(Iyi .* m22 + Izi .* m32), (GJ_Li .* m1 .* n1) + E4_Li.*(Iyi .* m2 .* n2 + Izi .* m3 .* n3);...
        (GJ_Li .* l1 .* n1) + E4_Li.*(Iyi .* l2 .* n2 + Izi .* l3 .* n3), (GJ_Li .* m1 .* n1) + E4_Li.*(Iyi .* m2 .* n2 + Izi .* m3 .* n3), (GJ_Li .* n12) + E4_Li.*(Iyi .* n22 + Izi .* n32)];
    
    k4 = [(-GJ_Li .* l12) + E2_Li.*(Iyi .* l22 + Izi .* l32), (-GJ_Li .* l1 .* m1) + E2_Li.*(Iyi .* l2 .* m2 + Izi .* l3 .* m3), (-GJ_Li .* l1 .* n1) + E2_Li.*(Iyi .* l2 .* n2 + Izi .* l3 .* n3);...
        (-GJ_Li .* l1 .* m1) + E2_Li.*(Iyi .* l2 .* m2 + Izi .* l3 .* m3), (-GJ_Li .* m12) + E2_Li.*(Iyi .* m22 + Izi .* m32), (-GJ_Li .* m1 .* n1) + E2_Li.*(Iyi .* m2 .* n2 + Izi .* m3 .* n3);...
        (-GJ_Li .* l1 .* n1) + E2_Li.*(Iyi .* l2 .* n2 + Izi .* l3 .* n3), (-GJ_Li .* m1 .* n1) + E2_Li.*(Iyi .* m2 .* n2 + Izi .* m3 .* n3), (-GJ_Li .* n12) + E2_Li.*(Iyi .* n22 + Izi .* n32)];
    
    k2T = permute(k2,[2 1 3]);
    
    k = [k1     k2     -k1      k2;...
        k2T    k3     -k2T     k4;...
        -k1    -k2      k1     -k2;...
        k2T    k4     -k2T     k3];
    
    dim = nBeams * DoF;
    globalStiffnessMatrix = zeros(dim,dim,1);
    ID = dim - (DoF * 2) + 1: dim;

    for j = nBeams:-1:1

        localStiffnessMatrix = k(:,:,i);
        
        if j == nBeams
            
            globalStiffnessMatrix(ID(7:12),ID(7:12),:) =  globalStiffnessMatrix(ID(7:12),ID(7:12),:) + localStiffnessMatrix(1:DoF,1:DoF,:);
            globalStiffnessMatrix(1:DoF,1:DoF,:) =  globalStiffnessMatrix(1:DoF,1:DoF,:) + localStiffnessMatrix(DoF+1:end,DoF+1:end,:);
            globalStiffnessMatrix(1:DoF,ID(7:12),:) =  globalStiffnessMatrix(1:DoF,ID(7:12),:) + localStiffnessMatrix(DoF+1:end,1:DoF,:);
            globalStiffnessMatrix(ID(7:12),1:DoF,:) =  globalStiffnessMatrix(ID(7:12),1:DoF,:) + localStiffnessMatrix(1:DoF,DoF+1:end,:);
        else
            globalStiffnessMatrix(ID,ID,:) =  globalStiffnessMatrix(ID,ID,:) + localStiffnessMatrix;
            ID = ID - DoF;
        end
    end
    
%     qf(13:end,i) = globalStiffnessMatrix(13:end,13:end,:)\loads(13:end,i);
    qf(7:end-6,i) = globalStiffnessMatrix(7:end-6,7:end-6,:)\loads(7:end-6,i);
end

displace(:,:,1) = qf(xForce,:);
displace(:,:,2) = qf(yForce,:);
displace(:,:,3) = qf(zForce,:);

rotation(:,:,1) = qf(xMoment,:);
rotation(:,:,2) = qf(yMoment,:);
rotation(:,:,3) = qf(zMoment,:);

displace = [displace(1:row,:,:), flipud(displace(row+1:end,:,:))];
rotation = [rotation(1:row,:,:), flipud(rotation(row+1:end,:,:))];

% Repeat for all leading/trailing edge panels
displace = [repmat(displace(1,:,:), first, 1); displace(2:end-1,:,:); repmat(displace(end,:,:), boxDim - last + 1, 1)];
rotation = [repmat(rotation(1,:,:), first, 1); rotation(2:end-1,:,:); repmat(rotation(end,:,:), boxDim - last + 1, 1)];

ua = 0.25 .* displace(1:end-1,1:end-1,:) + rotation(1:end-1,1:end-1,:) .* vec2(1:end-1,:,:)...
    + 0.25 .* displace(1:end-1,2:end,:) + rotation(1:end-1,2:end,:) .* vec1(1:end-1,:,:)...
    + 0.25 .* displace(2:end,1:end-1,:) + rotation(2:end,1:end-1,:) .* vec4(2:end,:,:)...
    + 0.25 .* displace(2:end,2:end,:) + rotation(2:end,2:end,:) .* vec3(2:end,:,:);

%% Check new wing
part.Points = part.Points + displace;
plotter(part);
%%

newCentre = wingCentre + ua;

plotter(part)
figure(2)
hold on
plot3(wingCentre(:,:,1),wingCentre(:,:,2),wingCentre(:,:,3),'r*')
plot3(newCentre(:,:,1),newCentre(:,:,2),newCentre(:,:,3),'r*')

end