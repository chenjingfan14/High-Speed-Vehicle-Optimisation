function [Cp,Mach,P] = prandtlmeyer(del,pm,Cp,Mach,P,flow,PrandtlMeyer)

Minf = flow.Minf;
Pinf = flow.Pinf;
gamma = flow.gamma;
vmax = (pi/2)*(((gamma+1)/(gamma-1))^0.5 - 1);

[row,col] = size(del);

dTheta = del(1:end-1,:) - del(2:end,:);

rowCon = any(pm,2);
rowArray = 1:row;
rowi = rowArray(rowCon);

nRows = size(PrandtlMeyer,1);

PM1 = PrandtlMeyer(:,1);
PM2 = PrandtlMeyer(:,2);

for i=rowi
    j = pm(i,:);
    colArray = (1:col)';
    colj = colArray(j);
    
    dThetai = dTheta(i-1,j)';
    
    id = i + (colj-1)*row;
    
    P1 = P(id-1);
    M1 = Mach(id-1);
    
    %% Should be shockwaves here
    % Although current vehicle gen should not make configs that have
    % contracting surfaces, may just be zero area panels causing negatives
    shock = dThetai < 0;
    
    P(id(shock)) = P1(shock);
    Mach(id(shock)) = M1(shock);
    
    %%
    con = M1 < 1 | shock;
    
    if all(con)
        continue
    end
    
    id(con) = [];
    P1(con) = [];
    M1(con) = [];
    dThetai(con) = [];
    
    dim = length(M1);
    
    %% Test Section
    
    % New and fast
    first = ones(dim,1);
    last = first*nRows;
    stopCon = false(dim,1);
    prevMiddle = zeros(dim,1);
    
    while any(~stopCon)
        middle = floor((first + last)/2);
        
        lThan = M1 < PM1(middle);
        last(lThan) = middle(lThan) - 1;
            
        gThan = M1 > PM1(middle);
        first(gThan) = middle(gThan) + 1;
            
        con = prevMiddle == middle;
        
        if any(con)
            stopCon(con) = true;
        end
        
        prevMiddle = middle;
        
    end
    
    rows = middle;
    
    % Half space search can be +-1 away from actual closest lookup value
    % Therefore interpolate between rows,rows-1 and rows+1,rows and average
    if any(rows == 1)
        vM1 = PM2(rows) + (M1 - PM1(rows)).*((PM2(rows+1) - PM2(rows))./(PM1(rows+1) - PM1(rows)));
    elseif any(rows == nRows)
        vM1 = PM2(rows-1) + (M1 - PM1(rows-1)).*((PM2(rows) - PM2(rows-1))./(PM1(rows) - PM1(rows-1)));
    else
        int1 = PM2(rows-1) + (M1 - PM1(rows-1)).*((PM2(rows) - PM2(rows-1))./(PM1(rows) - PM1(rows-1)));
        int2 = PM2(rows) + (M1 - PM1(rows)).*((PM2(rows+1) - PM2(rows))./(PM1(rows+1) - PM1(rows)));
        
        vM1 = (int1 + int2)/2;
    end

    vM2 = vM1 + dThetai;
    
    % New and fast
    first = ones(dim,1);
    last = first*nRows;
    stopCon = false(dim,1);
    prevMiddle = zeros(dim,1);
    
    while any(~stopCon)
        middle = floor((first + last)/2);
        
        lThan = vM2 < PM2(middle);
        last(lThan) = middle(lThan) - 1;
        
        gThan = vM2 > PM2(middle);
        first(gThan) = middle(gThan) + 1;
        
        con = prevMiddle == middle;
        
        if any(con)
            stopCon(con) = true;
        end
        
        prevMiddle = middle;
        
    end
    
    rows = middle;
    
    if any(rows == 1)
        M2 = PM1(rows) + (vM2 - PM2(rows)).*((PM1(rows+1) - PM1(rows))./(PM2(rows+1) - PM2(rows)));
    elseif any(rows == nRows)
        M2 = PM1(rows-1) + (vM2 - PM2(rows-1)).*((PM1(rows) - PM1(rows-1))./(PM2(rows) - PM2(rows-1)));
    else
        int1 = PM1(rows-1) + (vM2 - PM2(rows-1)).*((PM1(rows) - PM1(rows-1))./(PM2(rows) - PM2(rows-1)));
        int2 = PM1(rows) + (vM2 - PM2(rows)).*((PM1(rows+1) - PM1(rows))./(PM2(rows+1) - PM2(rows)));
        
        M2 = (int1 + int2)/2;
    end
    
    %%
    
    P2_P1 = ((1 + ((gamma-1)/2)*(M1.^2))./(1 + ((gamma-1)/2)*(M2.^2))).^(gamma/(gamma-1));
    P2 = P2_P1.*P1;
    
    Mach(id) = M2;
    P(id) = P2;
    
    Cp(id) = 2*((P2/Pinf)-1)/(gamma*(Minf^2));
    
    %% Alterations based on max turning angle
    % If turning angle greater than max, use different method to find panel
    % characteristics
    con = dThetai > vmax - vM1;
    
    %     % Base pressure
    %     Cp(i,con) = 1/(Minf^2);
    %     % Newtonian
    %     Cp(i,con) = 0;
    %     Mach(i,con) = Minf;
    %     P(i,con) = 0.5*Cp(i,con)*gamma*(Minf^2) + Pinf;
    
    % Vacuum
    Cp(i,con) = 0;
    Mach(i,con) = 0;
    P(i,con) = 0;
end

Cp = Cp(pm);
Mach = Mach(pm);
P = P(pm);