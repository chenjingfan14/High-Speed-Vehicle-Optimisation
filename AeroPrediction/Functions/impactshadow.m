function [impact,shadow] = impactshadow(pcy,pcz,area,prev,logicalFlow,yzUpBound,yzLoBound)

pcyArray = reshape(pcy,[],1);
pczArray = reshape(pcz,[],1);

% Check if panel lies within upper shield boundary
yUpCon = pcyArray <= yzUpBound(:,1)';
zUpCon = pczArray <= yzUpBound(:,2)';

% Check if panel lies within lower shield boundary
yLoCon = pcyArray <= yzLoBound(:,1)';
zLoCon = -pczArray <= -yzLoBound(:,2)';

% If between upper and lower shadow boundary, panel is shielded
% from freestream flow
shielded = any(yUpCon & zUpCon,2) & any(yLoCon & zLoCon,2);
shielded = reshape(shielded, size(logicalFlow));

nonZeroArea = area > 0;

% If panel normal faces flow, and has non-zero area, and is not
% shielded by previous panels, it is impacted directly by flow
impact = [prev; logicalFlow & ~shielded & nonZeroArea];
shadow = [prev; (~logicalFlow | shielded) & nonZeroArea];