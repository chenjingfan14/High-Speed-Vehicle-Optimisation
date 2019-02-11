function wing = wingstructure(wing,wingbox)

%% Wingbox x-Discretisation

sparMid = wingbox.Location;
sparThick = wingbox.SparThick;

di = wing.Dihedral;
wingPoints = wing.Points;
xNorm = wing.xNorm;
TE = wing.TrailingEdge;

sparBegin = (sparMid - sparThick/2);
sparEnd = (sparMid + sparThick/2);

numSpars = numel(sparBegin);

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
    
    [~,a] = min(abs(xNorm - sparBegin(i,:)));
    [~,b] = min(abs(xNorm - sparEnd(i,:)));
    
    if a == b
        
        b = a + 1;
    end
    
    if a == 1
        
        a = 2;
        b = 3;
    
    elseif b >= row
    
        b = row - 1;
        a = b - 1;
    end
    
    initSparPoints(:,:,1) = [sparBegin(i,:); sparEnd(i,:)] .* chords + xLE;
    xNorm([a b]) = [sparBegin(i,:); sparEnd(i,:)];
    
    for j = col:-1:1
        
        
        initSparPoints(:,j,2) = interp1(wingPoints(rows,j,1),wingPoints(rows,j,2),initSparPoints(:,j,1),'pchip');
        initSparPoints(:,j,3) = interp1(wingPoints(rows,j,1),wingPoints(rows,j,3),initSparPoints(:,j,1),'pchip');
    end
    
    wingPoints([a b],:,:) = initSparPoints;
    
    initSparPoints(:,1:half,[2 3]) = initSparPoints(:,1:half,[2 3]) - (skinThick(1:half) .* rotation);
    initSparPoints(:,half+1:end,[2 3]) = initSparPoints(:,half+1:end,[2 3]) + (skinThick(half+1:end) .* rotation);
    
    a = initSparPoints(:,1:half,:);
    b = fliplr(initSparPoints(:,half+1:end,:));
    
    for j = 3:-1:1
        
        % 4 points per spar
        sparPoints(:,:,j) = reshape([a(:,:,j) b(:,:,j)]', [], 4)';
    end
    
%     sparCell(i,:) = mat2cell(sparPoints, 4, ones(1,half), 3);
    sparCell{i} = sparPoints;
end

wing.Points = wingPoints;


sparEdges = [sparBegin(1); sparEnd(end)];
between = xNorm >= sparEdges(1) & xNorm <= sparEdges(2);

boxUpperOut = wingPoints(between,1:half,:);
boxLowerOut = wingPoints(between,half+1:end,:);

%% DIMENSIONALISE SKINTHICK

boxUpperIn = boxUpperOut;
boxLowerIn = boxLowerOut;

boxUpperIn(:,:,[2 3]) = boxUpperOut(:,:,[2 3]) - (skinThick(1:half) .* rotation);
boxLowerIn(:,:,[2 3]) = boxLowerOut(:,:,[2 3]) + (skinThick(half+1:end) .* rotation);

%% Proof plotting
figure
hold on
axis equal
plot3(wing.Points(:,:,1),wing.Points(:,:,2),wing.Points(:,:,3),'k')
plot3(wing.Points(:,:,1),wing.Points(:,:,2),wing.Points(:,:,3),'k.')

plot3(boxUpperOut(:,:,1),boxUpperOut(:,:,2),boxUpperOut(:,:,3),'r')
plot3(boxLowerOut(:,:,1),boxLowerOut(:,:,2),boxLowerOut(:,:,3),'r')

plot3(boxUpperIn(:,:,1),boxUpperIn(:,:,2),boxUpperIn(:,:,3),'r')
plot3(boxLowerIn(:,:,1),boxLowerIn(:,:,2),boxLowerIn(:,:,3),'r')

for i = 1:numel(sparCell)
    
    plot3(sparCell{i}(:,:,1),sparCell{i}(:,:,2),sparCell{i}(:,:,3),'r.')
    plot3(sparCell{i}(:,:,1),sparCell{i}(:,:,2),sparCell{i}(:,:,3),'r.')
end
hold off

for i = half:-1:1
    
    skin.Upper{i} = [boxUpperOut(:,i,:), boxUpperIn(:,i,:)];
    skin.Lower{i} = [boxLowerOut(:,i,:), boxLowerIn(:,i,:)];
    skin.UpperMean{i} = squeeze((boxUpperOut(:,i,:) + boxUpperIn(:,i,:))/2);
    skin.LowerMean{i} = squeeze((boxLowerOut(:,i,:) + boxLowerIn(:,i,:))/2);
end

wing.Box = between;
skin.Thickness = skinThick(1:half);
wing.Spar.Thickness = sparThick * chords(1:half);
wing.Skin = skin;
wing.Spar.Points = sparCell;

end