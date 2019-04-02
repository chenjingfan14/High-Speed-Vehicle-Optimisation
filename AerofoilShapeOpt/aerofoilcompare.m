function cost = aerofoilcompare(test,comparison)

% Looking for vertical comparison only
verticalError = abs(test(:,2) - comparison(:,2));

cost = mean(verticalError);% + max(verticalError);