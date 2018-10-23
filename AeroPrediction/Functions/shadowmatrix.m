function [yzUpBoundNew,yzLoBoundNew] = shadowmatrix(yi,zi,yzBound)

[row,col] = size(yi);
cols = 1:col;

[yMax,yMaxRows] = max(yi,[],1);
[zMax,zMaxRows] = max(zi,[],1);
[zMin,zMinRows] = min(zi,[],1);

yMaxID = yMaxRows + (cols-1)*row;
zMaxID = zMaxRows + (cols-1)*row;
zMinID = zMinRows + (cols-1)*row;

yziBound = [yMax,yi(zMaxID),yi(zMinID) ; zi(yMaxID),zMax,zMin]';
yziBound = unique(yziBound,'rows');

yzBound = sortrows([yzBound;yziBound],2,'descend');

[~,maxLoc] = max(yzBound(:,1));

yzUpBound = yzBound(1:maxLoc,:);
yzLoBound = yzBound(maxLoc:end,:);

yUpBound = yzUpBound(:,1);
zUpBound = yzUpBound(:,2);
yLoBound = yzLoBound(:,1);
zLoBound = yzLoBound(:,2);

ycon = yUpBound > yUpBound';
zcon = zUpBound > zUpBound';

upCon = ~any(ycon & zcon,1);

ycon = yLoBound > yLoBound';
zcon = -zLoBound > -zLoBound';

loCon = ~any(ycon & zcon,1);

yzUpBoundNew = yzUpBound(upCon,:);
yzLoBoundNew = yzLoBound(loCon,:);

yzBoundNew = [yzUpBoundNew;yzLoBoundNew];
% figure(3)
% hold on
% plot3(zeros(size(yzBoundNew,1)),yzBoundNew(:,1),yzBoundNew(:,2),'r*')
% hold off