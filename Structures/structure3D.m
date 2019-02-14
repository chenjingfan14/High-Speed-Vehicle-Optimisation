function structure(wing, part, wingP)
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

DoF = 6;

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
linTaperEq = 1 + spanTaper;

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

%% Chordwise beam stiffness properties
% Reshape points, put first node at end to complete circle
boxPoints = [boxPoints(transID); boxPoints(transID(1,:,:))];
boxNode1 = boxPoints(1:end-1,:,:);
boxNode2 = boxPoints(2:end,:,:);

% Create reference node
boxNode3 = (boxNode1 + boxNode2)/2;
r1 = boxNode2 - boxNode1;

chordVec = diff(boxPoints);

chordLength = (chordVec(:,:,1).^2 + chordVec(:,:,2).^2 + chordVec(:,:,3).^2).^0.5;

% No dihedral applied to inboard section, only needs to be altered in y dim
boxNode3(:,1,2) = boxNode3(:,1,2) + 0.5 * chordLength(:,1);
% node3(:,2:end,[2 3]) = node3(:,2:end,[2 3]) + 0.5 * permute(L(:,2:end),[1 3 2]) .* permute(repmat([cos(di); sin(di)], 1, nSecs-1), [3 2 1]);
boxNode3(:,2:end,[2 3]) = boxNode3(:,2:end,[2 3]) + 0.5 * chordLength(:,2:end) .* permute(repmat([cos(di); sin(di)], 1, nSecs-1), [3 2 1]);

r13 = boxNode3 - boxNode1;

r3(:,:,1) = r1(:,:,2).*r13(:,:,3) - r1(:,:,3).*r13(:,:,2);
r3(:,:,2) = r1(:,:,3).*r13(:,:,1) - r1(:,:,1).*r13(:,:,3);
r3(:,:,3) = r1(:,:,1).*r13(:,:,2) - r1(:,:,2).*r13(:,:,1);

r2(:,:,1) = r3(:,:,2).*r1(:,:,3) - r3(:,:,3).*r1(:,:,2);
r2(:,:,2) = r3(:,:,3).*r1(:,:,1) - r3(:,:,1).*r1(:,:,3);
r2(:,:,3) = r3(:,:,1).*r1(:,:,2) - r3(:,:,2).*r1(:,:,1);

r1Mag = (r1(:,:,1).^2 + r1(:,:,2).^2 + r1(:,:,3).^2).^0.5;
r2Mag = (r2(:,:,1).^2 + r2(:,:,2).^2 + r2(:,:,3).^2).^0.5;
r3Mag = (r3(:,:,1).^2 + r3(:,:,2).^2 + r3(:,:,3).^2).^0.5;

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

[r2, r3] = deal(zeros(size(r1)));

r3(:,:,1) = r1(:,:,2).*r13(:,:,3) - r1(:,:,3).*r13(:,:,2);
r3(:,:,2) = r1(:,:,3).*r13(:,:,1) - r1(:,:,1).*r13(:,:,3);
r3(:,:,3) = r1(:,:,1).*r13(:,:,2) - r1(:,:,2).*r13(:,:,1);

r2(:,:,1) = r3(:,:,2).*r1(:,:,3) - r3(:,:,3).*r1(:,:,2);
r2(:,:,2) = r3(:,:,3).*r1(:,:,1) - r3(:,:,1).*r1(:,:,3);
r2(:,:,3) = r3(:,:,1).*r1(:,:,2) - r3(:,:,2).*r1(:,:,1);

r1Mag = (r1(:,:,1).^2 + r1(:,:,2).^2 + r1(:,:,3).^2).^0.5;
r2Mag = (r2(:,:,1).^2 + r2(:,:,2).^2 + r2(:,:,3).^2).^0.5;
r3Mag = (r3(:,:,1).^2 + r3(:,:,2).^2 + r3(:,:,3).^2).^0.5;

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

% Initially set all radii as skin thickness
chordRad = repmat(skinThick,numel(index),1);
% Replace outer spar stations with spar thickness
chordRad(chordCon,:) = repmat(sparThick,numel(chordCon),1);
chordA = pi * chordRad.^2;

chordMass = chordA .* chordLength * rhoAl;

chordIy = 1/12 * chordMass .* (3 * chordRad .^2 + chordLength.^2);
chordIz = chordIy;

chordJ = chordIy + chordIz;

chordAE_L = chordA .* E ./ chordLength;
chordE12_L3 = 12 * E ./ chordLength.^3;
chordE6_L2 = 6 * E ./ chordLength.^2;
chordGJ_L = G * chordJ ./ chordLength;
chordE4_L = 4 * E ./ chordLength;
chordE2_L = 2 * E ./ chordLength;

spanRad = (sparThick(1:end-1) + sparThick(2:end))/2;
spanA = pi * spanRad.^2;

spanMass = spanA .* spanLength * rhoAl;

spanIy = 1/12 * spanMass .* (3 * spanRad .^2 + spanLength.^2);
spanIz = spanIy;

spanJ = spanIy + spanIz;

spanAE_L = linTaperEq .* spanA .* E ./ spanLength;
spanE12_L3 = 12 * E ./ spanLength.^3;
spanE6_L2 = 6 * E ./ spanLength.^2;
spanGJ_L = G * spanJ ./ spanLength;
spanE4_L = 4 * E ./ spanLength;
spanE2_L = 2 * E ./ spanLength;

%% Assemble load vector
xForce = 1:DoF:dim*2;
yForce = xForce + 1;
zForce = yForce + 1;
xMoment = zForce + 1;
yMoment = xMoment + 1;
zMoment = yMoment + 1;

[loads,qf] = deal(zeros(dim*2,1));

% Check order in which these are being applied to load vector
loads(xForce) = F(:,:,1);
loads(yForce) = F(:,:,2);
loads(zForce) = F(:,:,3);
loads(xMoment) = M(:,:,1);
loads(yMoment) = M(:,:,2);
loads(zMoment) = M(:,:,3);

% Complete circle so number of node 1s = number of total nodes
nodes = numel(boxNode1(:,:,1));
dim = nodes * DoF;
globalStiffnessMatrix = zeros(dim,dim,1);
ID = dim - (DoF * 2) + 1: dim;

% Matrix position arrays for local matrices
locBegin = 1:DoF;
locEnd = locBegin + DoF;

[chordBeams,~] = size(chordLength);

for i = nSecs:-1:1
    
    init = (i-1)*dim/nSecs + 1;
    
    kBox = elementalStiffnessMatrices(chordAE_L(:,i),chordE12_L3(:,i),chordE6_L2(:,i),...
        chordGJ_L(:,i),chordE4_L(:,i),chordE2_L(:,i),chordIy(:,i),chordIz(:,i),...
        boxlmn1(:,i,:),boxlmn2(:,i,:),boxlmn3(:,i,:));
    
    % Mid partitions will have two spar sets interefering, root and tip
    % only one
    if i ~= 1
        
        kSpar = elementalStiffnessMatrices(spanAE_L(:,i-1),spanE12_L3(:,i-1),spanE6_L2(:,i-1),...
            spanGJ_L(:,i-1),spanE4_L(:,i-1),spanE2_L(:,i-1),spanIy(:,i-1),spanIz(:,i-1),...
            sparlmn1(:,i-1,:),sparlmn2(:,i-1,:),sparlmn3(:,i-1,:));   
    end
    
    for j = chordBeams:-1:1
        
        localStiffnessMatrix = kBox(:,:,j);
        
        if i ~= 1 && index(j)
            
            kSparj = kSpar(:,:,index(j));
            
            %% Check spar connections here, potentially wrong nodes being selected
            % outboardNode = ID(1:DoF+1); % Never used?
            outboardNode = ID(DoF+1:end);
            inboardNode = outboardNode - dim/nSecs;
            
            globalStiffnessMatrix(inboardNode,inboardNode,:) = ...
                globalStiffnessMatrix(inboardNode,inboardNode,:) + kSparj(locBegin,locBegin,:);
            
            globalStiffnessMatrix(outboardNode,outboardNode,:) = ...
                globalStiffnessMatrix(outboardNode,outboardNode,:) + kSparj(locEnd,locEnd,:);
            
            globalStiffnessMatrix(outboardNode,inboardNode,:) = ...
                globalStiffnessMatrix(outboardNode,inboardNode,:) + kSparj(locEnd,locBegin,:);
            
            globalStiffnessMatrix(inboardNode,outboardNode,:) = ...
                globalStiffnessMatrix(inboardNode,outboardNode,:) + kSparj(locBegin,locEnd,:);
        end
            
        if j == chordBeams
            
            boxBegin = init:init + DoF - 1;
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
            
            ID = ID - DoF * 2;
        else
            globalStiffnessMatrix(ID,ID,:) =  globalStiffnessMatrix(ID,ID,:) + localStiffnessMatrix;
            ID = ID - DoF;
        end
    end
end

% First set of nodes at root are fixed, rest free
unconstrained = nodes/nSecs * DoF + 1 : dim;
qf(unconstrained) = globalStiffnessMatrix(unconstrained,unconstrained,:)\loads(unconstrained);

[reshape1,reshape2,~] = size(F);

displace(:,:,1) = reshape(qf(xForce),reshape1,reshape2);
displace(:,:,2) = reshape(qf(yForce),reshape1,reshape2);
displace(:,:,3) = reshape(qf(zForce),reshape1,reshape2);

rotation(:,:,1) = reshape(qf(xMoment),reshape1,reshape2);
rotation(:,:,2) = reshape(qf(yMoment),reshape1,reshape2);
rotation(:,:,3) = reshape(qf(zMoment),reshape1,reshape2);

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