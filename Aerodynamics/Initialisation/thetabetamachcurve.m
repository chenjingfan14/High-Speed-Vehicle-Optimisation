function beta = thetabetamachcurve(M,theta)
    
gamma = 1.4;

beta = (0:1:90)*pi/180;

theta = atan(2*cot(beta).*((M^2*(sin(beta).^2)-1)./(M^2*(gamma+cos(2*beta))+2)));
con = theta >= 0;

theta = theta(con);
beta = beta(con);

[~,iMax] = max(theta);

betaFine = beta(iMax-1):0.0001:beta(iMax+1);

thetaFine = atan(2*cot(betaFine).*((M^2*(sin(betaFine).^2)-1)./(M^2*(gamma+cos(2*betaFine))+2)));

[theta(i),~] = max(thetaFine);