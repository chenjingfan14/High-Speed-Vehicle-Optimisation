clear all
clc

l = 0.508;
l = 10;

xPanels = 50;
disc = 3;

x = (0:1/xPanels:1)';

switch disc
    case 1
        
    case 2
        x = 0.5*(1-cos(x*pi));
    case 3
        x = 1-cos(x*pi/2);
end

x = x*l;

yzPanels = 25;
theta = 0:pi/yzPanels:pi;

costFun = @completeaeroprediction;

load('nosebodyexperimental')
load('thetaBetaCurves.mat');

% experimentalData = experimentalData{5:6}2;

dim = length(experimentalData);
Aref = zeros(dim,1);

numericalCp = cell(dim,1);
numericalCN = cell(dim,1);
numericalCA = cell(dim,1);
numericalCl = cell(dim,1);
numericalCd = cell(dim,1);
numericalCm = cell(dim,1);
numericalImpact = cell(dim,1);
numericalShadow = cell(dim,1);

for Case = dim:-1:1
    
    con = x/l <= 0.45;
    con2 = x/l >= 0.45;
    
    xp = x(con);
    xp2 = x(con2);
    
    if any(Case == 5:6)
        con2 = 0.45 <= x/l & x/l <= 0.9;
        con3 = x/l >= 0.9;
        xp2 = x(con2);
        xp3 = x(con3);
        if xp3(1) ~= xp2(end)
            xp3 = [xp2(end);xp3];
        end
    end
    
    if xp2(1) ~= xp(end)
        xp2 = [xp(end);xp2];
    end
    
    switch Case
        case 1
            part = 0.16667*xp;
            part2 = repmat(0.075*l,size(xp2));
        case 2
            part = ((((1.3125^2)+(xp/l).*(0.9-(xp/l))).^0.5) - 1.3125)*l;
            part2 = repmat(0.075*l,size(xp2));
        case 3
            part = (0.1118*(xp/l).^0.5)*l;
            part2 = repmat(0.075*l,size(xp2));
        case 4
            part = ((((1.3125^2)+(xp/l).*(0.9-(xp/l))).^0.5) - 1.3125)*l;
            part2 = ((((14.41^2)-(((xp2/l)-0.45).^2)).^0.5) - 14.335)*l;
        case 5
            part = ((((1.3125^2)+(xp/l).*(0.9-(xp/l))).^0.5) - 1.3125)*l;
            part2 = 0.075*l;
            part3 = (0.075 - 0.105*((xp3/l)-0.9))*l;
        case 6
            part = ((((1.3125^2)+(xp/l).*(0.9-(xp/l))).^0.5) - 1.3125)*l;
            part2 = 0.075*l;
            part3 = (0.075 + 0.105*((xp3/l)-0.9))*l;
    end
    
    y = part.*sin(theta);
    z = part.*cos(theta);
    
    y2 = part2.*sin(theta);
    z2 = part2.*cos(theta);
    
    y2(1,:) = y(end,:);
    z2(1,:) = z(end,:);
    
    Properties.Conical = 1;
    
    if any(Case == 5:6)
        
        y3 = part3.*sin(theta);
        z3 = part3.*cos(theta);
        
        xp3(end+1) = xp3(end);
        y3(end+1,:) = 0;
        z3(end+1,:) = 0;
        
        y3(1,:) = y2(end,:);
        z3(1,:) = z2(end,:);
        
        config(3).Name = "aftbody";
        config(3).x = xp3;
        config(3).y = y3;
        config(3).z = z3;
        cellprop{3} = Properties;
        
    else
        
        xp2(end+1) = xp2(end);
        y2(end+1,:) = 0;
        z2(end+1,:) = 0;
    end
    
    config(1).Name = "forebody";
    config(1).x = xp;
    config(1).y = y;
    config(1).z = z;
    
    config(2).Name = "aftbody";
    config(2).x = xp2;
    config(2).y = y2;
    config(2).z = z2;
    
    cellprop{1} = Properties;
    cellprop{2} = Properties;
    
    config = xyztopoints(config);
    
    for i=1:length(cellprop)
        cellprop{i}.Points = config(i);
    end
    
    Aref(Case) = pi*max(part)^2;
    
    plotter(config) 
    
    configCell(Case,:) = cellprop;
    
end

% Can put in parallel
for Case = 1:dim
    
    MAC = 0;
    
    expData = experimentalData{Case};
    flow = flowparameters(expData(:,1));
    Mrange = [1:0.0001:10,10.1:0.1:100];
    PrandtlMeyer = prandtlmeyerlookup(Mrange,flow);
    
    [AoA,~,methodCN,methodCA,methodCl,methodCd,methodCm,impactMatrix,shadowMatrix] = feval(costFun,configCell(Case,:),flow,Aref(Case),MAC,thetaBetaM,maxThetaBetaM,PrandtlMeyer);
    
    numericalCp{Case} = methodCp;
    numericalCN{Case} = methodCN;
    numericalCA{Case} = methodCA;
    numericalCl{Case} = methodCl;
    numericalCd{Case} = methodCd;
    numericalCm{Case} = methodCm;
    numericalImpact{Case} = impactMatrix;
    numericalShadow{Case} = shadowMatrix;
end

for i = dim:-1:1
    expData = experimentalData{i};
    AoA = expData(:,1);
    
    CN = numericalCN{i};
    CA = numericalCA{i};
    Cl = numericalCl{i};
    Cd = numericalCd{i};
    Cm = numericalCm{i};
    impact = numericalImpact{i};
    shadow = numericalShadow{i};

    plotmethodsdata(AoA,CN,CA,Cl,Cd,Cm,impact,shadow,expData);
    
    CNdiff = abs(expData(:,2)' - CN);
    CAdiff = abs(expData(:,3)'/2 - CA);
    
    CNdiffbar(:,i) = mean(CNdiff,2);
    CAdiffbar(:,i) = mean(CAdiff,2);
    
    [nrow,~] = size(CN);
    numArray = (1:nrow)';
    
    CNsorted = sortrows([numArray,CNdiffbar(:,i)],2);
    CAsorted = sortrows([numArray,CAdiffbar(:,i)],2);
    
    CNsortedMethods(:,i) = CNsorted(:,1);
    CAsortedMethods(:,i) = CAsorted(:,1);
    
end

save('numericalNoseBody56.mat','AoA','numericalCN','numericalCA','numericalCl','numericalCd','numericalCm','numericalImpact','numericalShadow','CNsortedMethods','CAsortedMethods')
