function [displace, newWing, newPart] = structure3D(wing, part, wingP, H, L, LH)
% Initial structural program. See Low-fidelity aerostructural optimization
% of aircraft wings with a simplified wingbox model using OpenAeroStruct

% Material property
rhoAl = 2810;

if iscell(wing)
    
    wing = wing{:};
end

%% Initialise Skin
skin = wing.Skin;
nSecs = ceil(size(wingP,2)/2);
di = wing.Dihedral;

skinThick = skin.Thickness;

DoFelement = 6;

E = 70e9;
v = 0.29;
G = E / (2 * (1 + v));

%% Initiliase Spars (Only works for 4 points per spar)
sparNodes = wing.Spar.Nodes;

sparThick = wing.Spar.Thickness;

%% Total structural configuration properties (spanwise)

spanBeamVectors = diff(sparNodes,1,2);

spanLength = squeeze(spanBeamVectors(:,:,1).^2 + spanBeamVectors(:,:,2).^2 + spanBeamVectors(:,:,3).^2).^0.5;

spanTaper = (sparThick(1:end-1) - sparThick(2:end))./sparThick(2:end);

if isempty(spanTaper)
    
    linTaperEq = 1;
else
    linTaperEq = 1 + spanTaper;
end

%% Replace by creating box points earlier?
box = wing.Box;

first = find(box,1,'first');
last = find(box,1,'last');

[boxDim,~] = size(box);

%% CHECK: Points correspond to panel coordinates ahead of panel centre
% wingPoints = part.Points;
wingCentre = part.centre;
wingS = part.area;

boxPoints = skin.Nodes;
boxCentre = wingCentre(box,:,:);
boxP = wingP(box,:,:);
boxS = wingS(box,:);

structNode1 = boxPoints(:,1:end-1,:);
structNode2 = boxPoints(:,2:end,:);

structNode1 = [repmat(structNode1(1,:,:), first, 1); structNode1(2:end-1,:,:); repmat(structNode1(end,:,:), boxDim - last + 1, 1)];
structNode2 = [repmat(structNode2(1,:,:), first, 1); structNode2(2:end-1,:,:); repmat(structNode2(end,:,:), boxDim - last + 1, 1)];

force = wingP .* wingS;

vec1 = [wingCentre; structNode1(end,:,:)] - structNode1;
vec2 = [wingCentre; structNode2(end,:,:)] - structNode2;
vec3 = [structNode1(1,:,:); wingCentre] - structNode1;
vec4 = [structNode2(1,:,:); wingCentre] - structNode2;

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

crossLeft1 = crossmat(vec1sym, forceLeft1);
crossRight1 = crossmat(vec2sym, forceRight1);
crossLeft2 = crossmat(vec3sym, forceLeft2);
crossRight2 = crossmat(vec4sym, forceRight2);

M = 0.25 * crossLeft1 + 0.25 * crossRight1 +...
    0.25 * crossLeft2 + 0.25 * crossRight2 ;

%% Sum forces and moments at leading and trailing edge

F = [sum(F(1:first,:,:),1); F(first+1:last-1,:,:); sum(F(last:end,:,:),1)];
M = [sum(M(1:first,:,:),1); M(first+1:last-1,:,:); sum(M(last:end,:,:),1)];

Fs = H' * force(:);

% Rearrange wing spanwise wrap to wingbox chordwisel wrap
[row,col,~] = size(M);

F = reshape(Fs,row,col,3);

halfCol = ceil(col/2);

origID = zeros(row, col, 3);
origID(:) = 1:row*col*3;

transID = [origID(:, 1:halfCol, :); rot90(origID(:, halfCol + 1:end, :),2)];

F = F(transID);
M = M(transID);

M(:) = 0;

%% Chordwise beam stiffness properties
% Reshape points, put first node at end to complete circle
boxPoints = [boxPoints(transID); boxPoints(transID(1,:,:))];
boxNode1 = boxPoints(1:end-1,:,:);
boxNode2 = boxPoints(2:end,:,:);

% Create reference node
boxNode3 = (boxNode1 + boxNode2)/2;
r1 = boxNode2 - boxNode1;

chordVec = diff(boxPoints);

chordLength = magmat(chordVec);

% No dihedral applied to inboard section, only needs to be altered in y dim
boxNode3(:,1,2) = boxNode3(:,1,2) + 0.5 * chordLength(:,1);
% node3(:,2:end,[2 3]) = node3(:,2:end,[2 3]) + 0.5 * permute(L(:,2:end),[1 3 2]) .* permute(repmat([cos(di); sin(di)], 1, nSecs-1), [3 2 1]);
boxNode3(:,2:end,[2 3]) = boxNode3(:,2:end,[2 3]) + 0.5 * chordLength(:,2:end) .* permute(repmat([cos(di); sin(di)], 1, nSecs-1), [3 2 1]);

r13 = boxNode3 - boxNode1;

r3 = crossmat(r1, r13);
r2 = crossmat(r3, r1);

r1Mag = magmat(r1);
r2Mag = magmat(r2);
r3Mag = magmat(r3);

boxlmn1 = r1 ./ r1Mag;
boxlmn2 = r2 ./ r2Mag;
boxlmn3 = r3 ./ r3Mag;

%% Spanwise beam stiffness properties

sparNode1 = sparNodes(:,1:end-1,:);
sparNode2 = sparNodes(:,2:end,:);

% Create reference node upstream
sparNode3 = (sparNode1 + sparNode2)/2;
r1 = sparNode2 - sparNode1;

sparNode3(:,:,1) = sparNode3(:,:,1) + 0.5 * spanLength;

r13 = sparNode3 - sparNode1;

r3 = crossmat(r1, r13);
r2 = crossmat(r3, r1);

r1Mag = magmat(r1);
r2Mag = magmat(r2);
r3Mag = magmat(r3);

sparlmn1 = r1 ./ r1Mag;
sparlmn2 = r2 ./ r2Mag;
sparlmn3 = r3 ./ r3Mag;

%% Stiffness matrix properties

% Which box nodes are upper/lower spar rods attached to
% Able to use first chord/spanwise section as all sections will have same
% layout/node numbering

% Find way to avoid rounding here
a = permute(round(boxNode1(:,1,:),3),[1 3 2]);
b = permute(round(sparNodes(:,1,:),3),[1 3 2]);

% Use bottom spar nodes, as we are only looking for beam elements (ie one
% member per spar), whereas using both nodes would return two 
con = ismember(a,b(2:2:end,:,:),'rows');

% Only looking to use spars on the edge of the box here, inner spars will
% be incorporated within stiffness matrix loop
if sum(con)
    
    chordCon = [find(con,1,'first'); find(con,1,'last')];
else
    chordCon = con;
end

% Now where all spar nodes are contained within skin memebers
[~, index] = ismember(a,b,'rows');

%% Calculate elemental stiffness properties
% Initially set all radii as skin thickness
chordRad = repmat(skinThick,numel(index),1);
% Replace outer spar stations with spar thickness
chordRad(chordCon,:) = repmat(sparThick,numel(chordCon),1);
chordA = pi * chordRad.^2;

chordMass = chordA .* chordLength * rhoAl;

% chordIy = 1/12 * chordMass .* (3 * chordRad .^2 + chordLength.^2);
chordIy = pi/4 * chordRad.^4;
chordIz = chordIy;

chordJ = chordIy + chordIz;

chordAE_L = chordA .* E ./ chordLength;
chordE12_L3 = 12 * E ./ chordLength.^3;
chordE6_L2 = 6 * E ./ chordLength.^2;
chordGJ_L = G * chordJ ./ chordLength;
chordE4_L = 4 * E ./ chordLength;
chordE2_L = 2 * E ./ chordLength;

%% Update whether based on chord or not
% spanRad = (sparThick(1:end-1) + sparThick(2:end))/2;
spanRad = sparThick;
spanA = pi * spanRad.^2;

spanMass = spanA .* spanLength * rhoAl;
% spanIy = 1/12 * spanMass .* (3 * spanRad .^2 + spanLength.^2);

spanIy = pi/4 * spanRad.^4;
spanIz = spanIy;

spanJ = spanIy + spanIz;

spanAE_L = linTaperEq .* spanA .* E ./ spanLength;
spanE12_L3 = 12 * E ./ spanLength.^3;
spanE6_L2 = 6 * E ./ spanLength.^2;
spanGJ_L = G * spanJ ./ spanLength;
spanE4_L = 4 * E ./ spanLength;
spanE2_L = 2 * E ./ spanLength;

[chordBeams,~] = size(chordLength);

array = 1:chordBeams;

array = array(index > 0);

%% Mass properties

lumpedSparMass = [spanMass(:,1)/2 ...
    (spanMass(:,1:end-1) + spanMass(:,2:end))/2 ...
    spanMass(:,end)/2];

for i = chordBeams:-1:1
    
    if i == 1
        
        mass(i,:) = (chordMass(1,:) + chordMass(end,:))/2; 
    else
        mass(i,:) = (chordMass(i-1,:) + chordMass(i,:))/2;
    end
        
    if index(i)
    
        mass(i,:) = mass(i,:) + lumpedSparMass(index(i),:);
    end
end

%% Assemble load vector
nodes = numel(boxNode1(:,:,1));
DoFtotal = nodes * DoFelement;

xForce = 1:DoFelement:DoFtotal;
yForce = xForce + 1;
zForce = yForce + 1;
xMoment = zForce + 1;
yMoment = xMoment + 1;
zMoment = yMoment + 1;

[loads,qf] = deal(zeros(DoFtotal,1));

% Check order in which these are being applied to load vector
loads(xForce) = F(:,:,1);
loads(yForce) = F(:,:,2);
loads(zForce) = F(:,:,3) - mass * 9.81;
loads(xMoment) = M(:,:,1);
loads(yMoment) = M(:,:,2);
loads(zMoment) = M(:,:,3);

% Complete circle so number of node 1s = number of total nodes
globalStiffnessMatrix = zeros(DoFtotal,DoFtotal,1);
ID = DoFtotal - (DoFelement * 2) + 1: DoFtotal;

% Matrix position arrays for local matrices
locBegin = 1:DoFelement;
locEnd = locBegin + DoFelement;

for i = nSecs:-1:1
    
    init = (i-1)*DoFtotal/nSecs + 1;
    
%     kBox = elementalStiffnessMatrices(chordAE_L(:,i),chordE12_L3(:,i),chordE6_L2(:,i),...
%         chordGJ_L(:,i),chordE4_L(:,i),chordE2_L(:,i),chordIy(:,i),chordIz(:,i),...
%         boxlmn1(:,i,:),boxlmn2(:,i,:),boxlmn3(:,i,:));
    
    kBox = elementalStiffnessMatrices(chordAE_L(:,i),chordE12_L3(:,i),chordE6_L2(:,i),...
        chordGJ_L(:,i),chordE4_L(:,i),chordE2_L(:,i),chordIy,chordIz,...
        boxlmn1(:,i,:),boxlmn2(:,i,:),boxlmn3(:,i,:));
    
    % Test to check if all local stiffness matrices are symmetric and
    % positive definite
    flag = checkstiffnessmatrix(kBox);
    
    % Mid partitions will have two spar sets interefering, root and tip
    % only one
    if i ~= 1
        
        kSpar = elementalStiffnessMatrices(spanAE_L(:,i-1),spanE12_L3(:,i-1),spanE6_L2(:,i-1),...
            spanGJ_L(:,i-1),spanE4_L(:,i-1),spanE2_L(:,i-1),spanIy,spanIz,...
            sparlmn1(:,i-1,:),sparlmn2(:,i-1,:),sparlmn3(:,i-1,:));   
    end
    
    flag = checkstiffnessmatrix(kSpar);
    
    for j = chordBeams:-1:1
        
        localStiffnessMatrix = kBox(:,:,j);
        
        if i ~= 1 && index(j)
            
            kSparj = kSpar(:,:,index(j));
            
            %% Check spar connections here, potentially wrong nodes being selected
            % outboardNode = ID(1:DoF+1); % Never used?
            outboardNode = ID(DoFelement+1:end);
            inboardNode = outboardNode - DoFtotal/nSecs;
            
%             globalStiffnessMatrix(inboardNode,inboardNode,:) = ...
%                 globalStiffnessMatrix(inboardNode,inboardNode,:) + kSparj(locBegin,locBegin,:);
%             
%             globalStiffnessMatrix(outboardNode,outboardNode,:) = ...
%                 globalStiffnessMatrix(outboardNode,outboardNode,:) + kSparj(locEnd,locEnd,:);
            
            globalStiffnessMatrix(outboardNode,inboardNode,:) = ...
                globalStiffnessMatrix(outboardNode,inboardNode,:) + kSparj(locEnd,locBegin,:);
            
            globalStiffnessMatrix(inboardNode,outboardNode,:) = ...
                globalStiffnessMatrix(inboardNode,outboardNode,:) + kSparj(locBegin,locEnd,:);
        end
            
        if j == chordBeams
            
            boxBegin = init:init + DoFelement - 1;
            boxEnd = ID(7:12);
            
            globalStiffnessMatrix(boxEnd,boxEnd,:) = ...
                globalStiffnessMatrix(boxEnd,boxEnd,:) + localStiffnessMatrix(locBegin,locBegin,:);
            
            globalStiffnessMatrix(boxBegin,boxBegin,:) = ...
                globalStiffnessMatrix(boxBegin,boxBegin,:) + localStiffnessMatrix(locEnd,locEnd,:);
            
            globalStiffnessMatrix(boxBegin,boxEnd,:) = ...
                globalStiffnessMatrix(boxBegin,boxEnd,:) + localStiffnessMatrix(locEnd,locBegin,:);
            
            globalStiffnessMatrix(boxEnd,boxBegin,:) = ...
                globalStiffnessMatrix(boxEnd,boxBegin,:) + localStiffnessMatrix(locBegin,locEnd,:);
            
        elseif j == 1
            
            ID = ID - DoFelement * 2;
        else
            globalStiffnessMatrix(ID,ID,:) =  globalStiffnessMatrix(ID,ID,:) + localStiffnessMatrix;
            ID = ID - DoFelement;
        end
    end
end

%% Solve system of linear equations
% First set of nodes at root are fixed, rest free
unconstrained = nodes/nSecs * DoFelement + 1 : DoFtotal;

opts.SYM = true;

% Check if stiffness matrix is positive definite
% [~,p] = chol(globalStiffnessMatrix);

qf(unconstrained) = linsolve(globalStiffnessMatrix(unconstrained,unconstrained,:),loads(unconstrained),opts);

% globalStiffnessMatrix(unconstrained,unconstrained,:) = globalStiffnessMatrix(unconstrained,unconstrained,:)*10;
% qf = globalStiffnessMatrix\loads;

[reshape1,reshape2,~] = size(F);

j = 1:3;
for i = 1:6:length(qf)
    
    displace(j,:) = qf(i:i+2);
    j = j + 3;
end

[dim1,dim2,~] = size(skin.Nodes);

displace3D(:,:,1) = reshape(qf(xForce),dim1,dim2);
displace3D(:,:,2) = reshape(qf(yForce),dim1,dim2);
displace3D(:,:,3) = reshape(qf(zForce),dim1,dim2);

skinNodes3D = skin.Nodes + displace3D;
skinNodes1D = skin.Nodes(:) + displace;
    
% newCentre = H * skinNodes1D;
newPoints = LH * skinNodes1D;

% newCentre = reshape(newCentre,size(part.centre));
newPoints = reshape(newPoints,size(part.Points));

plotter(part)
figure(2)
hold on
plot3(newPoints(:,:,1),newPoints(:,:,2),newPoints(:,:,3),'k*')
% plot3(wingCentre(:,:,1),wingCentre(:,:,2),wingCentre(:,:,3),'r*')
% plot3(newCentre(:,:,1),newCentre(:,:,2),newCentre(:,:,3),'k*')

plotter(newPoints);

newPart = normals(newPoints);

newWing = wing;
newWing.Skin.Nodes = skinNodes3D;

%% Transferring skin nodes to make corresponding spar nodes
% Way to just flag skin node as spar instead of making separate structure?
half = floor(size(skinNodes3D,2)/2);
skinToSpar = [skinNodes3D(:,1:half,:); rot90(skinNodes3D(:,half+1:end,:),2)];

index(index == 0) = [];

newWing.Spar.Nodes(index,:,:) = skinToSpar(array,:,:);

close all

end