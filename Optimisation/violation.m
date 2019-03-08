function penalty = violation(constrain,minVal,maxVal)

penalty = zeros(size(constrain));

eta = 10;

tooLow = constrain < minVal;
tooBig = constrain > maxVal;

p1 = (minVal - constrain)./(maxVal - minVal);
p2 = (constrain - maxVal)./(maxVal - minVal);

penalty(tooLow) = p1(tooLow);
penalty(tooBig) = p2(tooBig);

penalty = eta*(sum(penalty,2)).^2;