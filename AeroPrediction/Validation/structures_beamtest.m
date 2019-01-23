% function structures_beamtest(beam)
% Initial structural program. See Low-fidelity aerostructural optimization
% of aircraft wings with a simplified wingbox model using OpenAeroStruct

%% Initialise Skin
skin = beam.Skin;

skinThick = skin.Thickness;
skinUp = skin.Upper.xyz;
skinLo = skin.Lower.xyz;

xSkinUpOut = skin.Upper.x(:,1,:);
xSkinLoOut = skin.Lower.x(:,1,:);
xSkinUpIn = skin.Upper.x(:,2,:);
xSkinLoIn = skin.Lower.x(:,2,:);

ySkinUpOut = skin.Upper.y(:,1,:);
ySkinLoOut = skin.Lower.y(:,1,:);
ySkinUpIn = skin.Upper.y(:,2,:);
ySkinLoIn = skin.Lower.y(:,2,:);

zSkinUpOut = skin.Upper.z(:,1,:);
zSkinLoOut = skin.Lower.z(:,1,:);
zSkinUpIn = skin.Upper.z(:,2,:);
zSkinLoIn = skin.Lower.z(:,2,:);


%% Calculate skin mid-line properties
% Find middle points of skins, use inner second and second last twice 
% instead of first and last to account for spar thickness

xMidUpper = (xSkinUpOut + xSkinUpIn)/2;
yMidUpper = (ySkinUpOut + ySkinUpIn)/2;
zMidUpper = (zSkinUpOut + zSkinUpIn)/2;

xMidLower = (xSkinLoOut + xSkinLoIn)/2;
yMidLower = (ySkinLoOut + ySkinLoIn)/2;
zMidLower = (zSkinLoOut + zSkinLoIn)/2;

ax = xMidUpper(1:end-1,:);
bx = xMidUpper(2:end,:);
cx = xMidLower(2:end,:);
dx = xMidLower(1:end-1,:);

ay = yMidUpper(1:end-1,:);
by = yMidUpper(2:end,:);
cy = yMidLower(2:end,:);
dy = yMidLower(1:end-1,:);

az = zMidUpper(1:end-1,:);
bz = zMidUpper(2:end,:);
cz = zMidLower(2:end,:);
dz = zMidLower(1:end-1,:);

skinUpperL = ((bx - ax).^2 + (by - ay).^2 + (bz - az).^2) .^0.5;
skinLowerL = ((cx - dx).^2 + (cy - dy).^2 + (cz - dz).^2) .^0.5;
skinIntL_t = sum([skinUpperL; skinLowerL]./skinThick,1);

%% Calculate skin parallelogram areas and trapezoidal box areas
% Areas based on xyz vectors so any triangle can be input
% Triangle one
skinUpperArea = skinUpperL .* skinThick;

height1 = sqrt((dx - ax).^2 + (dy - ay).^2 + (dz - az).^2);
A1 = 0.5 .* skinUpperL .* height1;

% Triangle two
skinLowerArea = skinLowerL .* skinThick;

height2 = sqrt((cx - bx).^2 + (cy - by).^2 + (cz - bz).^2);
A2 = 0.5 .* skinLowerL .* height2;

Ae = sum(A1 + A2,1);

%% Calculate skin centroid/MOI contributions

skinUpperCentroid = panelcentroid(skinUp);
skinLowerCentroid = panelcentroid(skinLo);

xSkinUpperCentroid = skinUpperCentroid(:,1:3:end);
xSkinLowerCentroid = skinLowerCentroid(:,1:3:end);

ySkinUpperCentroid = skinUpperCentroid(:,2:3:end);
ySkinLowerCentroid = skinLowerCentroid(:,2:3:end);

zSkinUpperCentroid = skinUpperCentroid(:,3:3:end);
zSkinLowerCentroid = skinLowerCentroid(:,3:3:end);

totSkinArea = sum([skinUpperArea; skinLowerArea],1);

skinUpperAx = skinUpperArea .* xSkinUpperCentroid;
skinLowerAx = skinLowerArea .* xSkinLowerCentroid;

skinUpperAy = skinUpperArea .* ySkinUpperCentroid;
skinLowerAy = skinLowerArea .* ySkinLowerCentroid;

skinUpperAz = skinUpperArea .* zSkinUpperCentroid;
skinLowerAz = skinLowerArea .* zSkinLowerCentroid;

%% Calculate skin moments of inertia
xDiffUpper = diff(xSkinUpOut);
zDiffUpper = diff(zSkinUpIn);

aUpper = zDiffUpper./xDiffUpper;
bUpper = (zDiffUpper + skinThick)/2;
cUpper = xDiffUpper;

xDiffLower = diff(xSkinLoOut);
zDiffLower = diff(zSkinLoIn);

aLower = -zDiffLower./xDiffLower;
bLower = (-zDiffLower + skinThick)/2;
cLower = zDiffLower;

% Momemnts of inertia about centroid
IpxSkinUpper = 2 * ((aUpper.^3 .* cUpper.^4)/12 + (aUpper.^2 .* bUpper .* cUpper.^3)/3 +...
    (aUpper .* bUpper.^2 .* cUpper.^2)/2 + (bUpper.^3 .* cUpper)/3);

IpxSkinLower = 2 * ((aLower.^3 .* cLower.^4)/12 + (aLower.^2 .* bLower .* cLower.^3)/3 +...
    (aLower .* bLower.^2 .* cLower.^2)/2 + (bLower.^3 .* cLower)/3);

IpzSkinUpper = skinThick .* skinUpperL;

IpzSkinLower = skinThick .* skinLowerL;

%% Initiliase Spars (Only works for 4 points per spar)
spars = beam.Spars.Points;

nSpars = numel(spars);
nSecs = beam.Partitions + 1;
sparPerSec = nSpars/nSecs;

sparThick = beam.Spars.Thickness;

% sparFront = spars.Front.xyz;
% sparBack = spars.Back.xyz;
% 
% sparMid = (sparFront + sparBack)/2;
% 
% dim = size(spars);
% sparMid = cell(dim);



for i = nSpars:-1:1

    spar = spars{i};
    
    sparFront(:,:,i) = spar(:,1:3);
    sparBack(:,:,i) = spar(:,end-2:end);
    sparMid(:,:,i) = (sparFront(:,:,i) + sparBack(:,:,i))/2;
    
    sparCentroid(i,:) = panelcentroid(spar);
end

%% Calculate spar mid-line properties
sparMidDiff = diff(sparMid);

sparMidHeight = (sparMidDiff(:,1,:).^2 + sparMidDiff(:,2,:).^2 + sparMidDiff(:,3,:).^2).^0.5;
sparMidHeight = reshape(sparMidHeight, sparPerSec, nSecs);

sparIntL_t = sum(sparMidHeight ./ sparThick,1);

sparFrontDiff = diff(sparFront);
sparBackDiff = diff(sparBack);

frontHeight = (sparFrontDiff(:,1,:).^2 + sparFrontDiff(:,2,:).^2 + sparFrontDiff(:,3,:).^2) .^0.5;
backHeight = (sparBackDiff(:,1,:).^2 + sparBackDiff(:,2,:).^2 + sparBackDiff(:,3,:).^2) .^0.5;

sparHeight = [frontHeight, backHeight];

aSpar = min(sparHeight,[],2);
bSpar = max(sparHeight,[],2);

aSpar = reshape(aSpar, sparPerSec, nSecs);
bSpar = reshape(bSpar, sparPerSec, nSecs);
%% Calculate spar centroid/MOI contributions
% For centroid we only need outer spar properties
sparArea = sparThick .* (aSpar + bSpar)/2;

outerArea = sparArea;
totOuterSparArea = sum(outerArea,1);

dim = size(outerArea);

xSparCentroid = reshape(sparCentroid(:,1:3:end),dim);
ySparCentroid = reshape(sparCentroid(:,2:3:end),dim);
zSparCentroid = reshape(sparCentroid(:,3:3:end),dim);

xOuterSparCentroid = xSparCentroid;
yOuterSparCentroid = ySparCentroid;
zOuterSparCentroid = zSparCentroid;

sparAx = outerArea .* xOuterSparCentroid;
sparAy = outerArea .* yOuterSparCentroid;
sparAz = outerArea .* zOuterSparCentroid;

% Momemnts of inertia about centroid
% Check directions?
% About xc axis
IpxSpar = (sparThick .* (aSpar + bSpar) .* (aSpar.^2 + bSpar.^2))/48;
% About zc axis
IpzSpar = (sparThick.^3 .* (aSpar.^2 + (4 * aSpar .* bSpar) + bSpar.^2)) ./...
    (36 .* (aSpar + bSpar));

%% Total structural configuration properties

intL_t = sparIntL_t + skinIntL_t;

J = (4 * Ae.^2) ./ intL_t;

shearCentre = [mean(xOuterSparCentroid);...
    mean(yOuterSparCentroid);...
    mean(zOuterSparCentroid)]';

% FIX THIS
xc = 1./(totOuterSparArea + totSkinArea) .* (sum([skinUpperAx; skinLowerAx],1) + sum(sparAx,1));
yc = 1./(totOuterSparArea + totSkinArea) .* (sum([skinUpperAy; skinLowerAy],1) + sum(sparAy,1));
zc = 1./(totOuterSparArea + totSkinArea) .* (sum([skinUpperAz; skinLowerAz],1) + sum(sparAz,1));

dxySkinUpper = ((xSkinUpperCentroid - xc).^2 + (ySkinUpperCentroid - yc).^2).^0.5;
dyzSkinUpper = ((ySkinUpperCentroid - yc).^2 + (zSkinUpperCentroid - zc).^2).^0.5;
dxzSkinUpper = ((xSkinUpperCentroid - xc).^2 + (zSkinUpperCentroid - zc).^2).^0.5;

dxySkinLower = ((xSkinLowerCentroid - xc).^2 + (ySkinLowerCentroid - yc).^2).^0.5;
dyzSkinLower = ((ySkinLowerCentroid - yc).^2 + (zSkinLowerCentroid - zc).^2).^0.5;
dxzSkinLower = ((xSkinLowerCentroid - xc).^2 + (zSkinLowerCentroid - zc).^2).^0.5;

dxySpar = ((xSparCentroid - xc).^2 + (ySparCentroid - yc).^2).^0.5;
dyzSpar = ((ySparCentroid - yc).^2 + (zSparCentroid - zc).^2).^0.5;
dxzSpar = ((xSparCentroid - xc).^2 + (zSparCentroid - zc).^2).^0.5;

Ipx = sum(IpxSkinUpper + skinUpperArea .* dyzSkinUpper.^2) + ...
    sum(IpxSkinLower + skinLowerArea .* dyzSkinLower.^2) + ...
    sum(IpxSpar + sparArea .* dyzSpar.^2);

% Ipy = sum((IpySkinUpper + skinUpperArea .* dxzSkinUpper.^2) + ...
%     (IpySkinLower + skinLowerArea .* dxzSkinLower.^2) + ...
%     (IpySpar + sparArea .* dxzSpar.^2));

Ipz = sum(IpzSkinUpper + skinUpperArea .* dxySkinUpper.^2) + ...
    sum(IpzSkinLower + skinLowerArea .* dxySkinLower.^2) + ...
    sum(IpzSpar + sparArea .* dxySpar.^2);

%% Plot
figure
hold on
axis equal
plot3(skinUp(:,1:3:end),skinUp(:,2:3:end),skinUp(:,3:3:end),'r')
% plot3(xMidUpper,yMidUpper,zMidUpper,'r--')

plot3(skinLo(:,1:3:end),skinLo(:,2:3:end),skinLo(:,3:3:end),'r')
% plot3(xMidLower,yMidLower,zMidLower,'r--')

plot3(skinUpperCentroid(:,1:3:end),skinUpperCentroid(:,2:3:end),skinUpperCentroid(:,3:3:end),'k.')
plot3(skinLowerCentroid(:,1:3:end),skinLowerCentroid(:,2:3:end),skinLowerCentroid(:,3:3:end),'k.')

plot3(sparCentroid(:,1:3:end),sparCentroid(:,2:3:end),sparCentroid(:,3:3:end),'k.')

for i = 1:nSpars
    % 1:3:end here ensures 3D matrices plot ALL spars
    spar = spars{i};
    sparMid = mean(spar);
    
    plot3(spar(:,1:3:end)',spar(:,2:3:end)',spar(:,3:3:end)','r')
%     plot3(sparMid(:,1:3:end)',sparMid(:,2:3:end)',sparMid(:,3:3:end)','r--')
end
plot3(xc,yc,zc,'bx')

hold off
