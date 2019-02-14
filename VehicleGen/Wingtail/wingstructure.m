function wing = wingstructure(wing,wingbox)

%% Wingbox x-Discretisation

sparMid = wingbox.Location;
sparThick = wingbox.SparThick;

di = wing.Dihedral;
wingPoints = wing.Points;
xNorm = wing.xNorm;
TE = wing.TrailingEdge;

numSpars = numel(sparMid);

% xNorm = sort([xNorm; sparBegin; sparEnd]);
chordVec = diff(wingPoints([1 end],:,:));
chords = (chordVec(:,:,1).^2 + chordVec(:,:,2).^2 + chordVec(:,:,3).^2).^0.5;
xLE = wingPoints(1,:,1);

skinThick = wingbox.SkinThick * chords;

[row, col, ~] = size(wingPoints);

% If wing has non-zero trailing edge, do not allow spars to interfere
if TE
    
    row = row - 1;
end

rows = 1:row - 1;

half = col/2;

% Create y, z rotation matrix, centreline chord should
% have zero degree rotation for symmetry
rotation = [0, 1; repmat([-sin(di) cos(di)], half-1, 1)];
rotation = permute(rotation,[3,1,2]);

for i = numSpars:-1:1
    
    [~,a] = min(abs(xNorm - sparMid(i,:)));
    
    initSparPoints(:,:,1) = sparMid(i,:) .* chords + xLE;
    xNorm(a) = sparMid(i);
    
    for j = col:-1:1
        
        initSparPoints(:,j,2) = interp1(wingPoints(rows,j,1),wingPoints(rows,j,2),initSparPoints(:,j,1),'pchip');
        initSparPoints(:,j,3) = interp1(wingPoints(rows,j,1),wingPoints(rows,j,3),initSparPoints(:,j,1),'pchip');
    end
    
    wingPoints(a,:,:) = initSparPoints;
    
    initSparPoints(:,1:half,[2 3]) = initSparPoints(:,1:half,[2 3]) - (skinThick(1:half)/2 .* rotation);
    initSparPoints(:,half+1:end,[2 3]) = initSparPoints(:,half+1:end,[2 3]) + (skinThick(half+1:end)/2 .* rotation);
    
    topSpar = initSparPoints(:,1:half,:);
    botSpar = fliplr(initSparPoints(:,half+1:end,:));
    
    for j = 3:-1:1
        % 2 points per spar
        sparPoints(:,:,j) = reshape([topSpar(:,:,j) botSpar(:,:,j)]', [], 2)';
    end
    
    sparNodes(i * 2 - 1,:,:) = mean(topSpar,1);
    sparNodes(i * 2,:,:) = mean(botSpar,1);
    sparCell{i} = sparPoints;
end

wing.Points = wingPoints;

between = xNorm >= sparMid(1) & xNorm <= sparMid(end);

boxOut = wingPoints(between,:,:);
boxIn = boxOut;

boxIn(:,1:half,[2 3]) = boxOut(:,1:half,[2 3]) - (skinThick(1:half) .* rotation);
boxIn(:,half+1:end,[2 3]) = boxOut(:,half+1:end,[2 3]) + (skinThick(half+1:end) .* rotation);

% for i = half:-1:1
% 
%     skin.Nodes{i} = squeeze((boxOut(:,i,:) + boxIn(:,i,:))/2);
% end

skinNodes = (boxOut + boxIn)/2;

skin.Nodes = skinNodes;
skin.Thickness = skinThick(1:half) .* chords(1:half);
wing.Spar.Thickness = sparThick * chords(1:half);
wing.Box = between;
wing.Skin = skin;
wing.Spar.Points = sparCell;
wing.Spar.Nodes = sparNodes;

%% Proof plotting
figure
hold on
axis equal
plot3(wing.Points(:,:,1),wing.Points(:,:,2),wing.Points(:,:,3),'k')

plot3(boxOut(:,:,1),boxOut(:,:,2),boxOut(:,:,3),'r')
plot3(boxIn(:,:,1),boxIn(:,:,2),boxIn(:,:,3),'r')
% plot3(skin.Nodes(:,:,1),skin.Nodes(:,:,2),skin.Nodes(:,:,3),'r.')

for i = 1:numel(sparCell)
    
    plot3(sparCell{i}(:,:,1),sparCell{i}(:,:,2),sparCell{i}(:,:,3),'r.')
    plot3(sparCell{i}(:,:,1),sparCell{i}(:,:,2),sparCell{i}(:,:,3),'r.')
end
hold off

end