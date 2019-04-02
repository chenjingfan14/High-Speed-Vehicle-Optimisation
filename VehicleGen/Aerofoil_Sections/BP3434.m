function foilData = BP3434(variables,varArray,n,nPop,disc,lDisc)

if nargin >= 6
    
    Pxu = disc;
    Pxl = lDisc;
    
elseif nargin >= 5
    
    Pxu = (0:1/(disc - 1):1)';
    Pxl = Pxu;
else
    Pxu = (0:1/100:1)';
    Pxl = Pxu;
end

gleArray = varArray == "gle";
b0Array = varArray == "b0";
b2Array = varArray == "b2";
xcArray = varArray == "xc";
ycArray = varArray == "yc";
b17Array = varArray == "b17";
zteArray = varArray == "zte";
ateArray = varArray == "ate";
rleArray = varArray == "rle";
b8Array = varArray == "b8";
xtArray = varArray == "xt";
ytArray = varArray == "yt";
b15Array = varArray == "b15";
dzteArray = varArray == "dzte";
bteArray = varArray == "bte";
b1Array = varArray == "b1";

foilData = cell(nPop,n + 1);

flag = false(nPop,1);

for i = 1:nPop
    
    for j = 1:n + 1
        
        gle = variables(i, gleArray);
        b0 = variables(i, b0Array);
        b2 = variables(i, b2Array);
        xc = variables(i, xcArray);
        yc = variables(i, ycArray);
        b17 = variables(i, b17Array);
        zte = variables(i, zteArray);
        ate = variables(i, ateArray);
        rle = variables(i, rleArray);
        b8 = variables(i, b8Array);
        xt = variables(i, xtArray);
        yt = variables(i, ytArray);
        b15 = variables(i, b15Array);
        dzte = variables(i, dzteArray);
        bte = variables(i, bteArray);
        b1 = variables(i, b1Array);
        
        % Discretisation will depend on xt, xc locations as these are where curves
        % merge
        
        nDisc = 100;
        
        uc = round(xc * nDisc);
        ut = round(xt * nDisc);
        
        ulc = (0 : 1/(uc - 1) : 1)';
        utc = (0 : 1/(nDisc - uc + 1) : 1)';
        ult = (0 : 1/(ut - 1) : 1)';
        utt = (0 : 1/(nDisc - ut + 1) : 1)';
        
        %% BP3434 Paramters
        
        cgle = cot(gle);
        cate = cot(ate);
        
        %% 3 degree Bezier curve
        % Leading thickness curve
        
        upperLimit = min(yt, (-2*rle*xt/3).^0.5);
        con = b8 >= upperLimit;
        
        b8(b8 <= 0) = 0.01;
        b8(con) = upperLimit(con) - 0.01;
        
        xlt0 = 0;
        xlt1 = 0;
        xlt2 = (-3 * b8.^2)/(2 * rle);
        xlt3 = xt;
        
        ylt0 = 0;
        ylt1 = b8;
        ylt2 = yt;
        ylt3 = yt;
        
        xlt = [xlt0 xlt1 xlt2 xlt3];
        ylt = [ylt0 ylt1 ylt2 ylt3];
        
        % Leading camber curve
        xlc0 = 0;
        xlc1 = b0;
        xlc2 = b2;
        xlc3 = xc;
        
        ylc0 = 0;
        ylc1 = b0 * tan(gle);
        ylc2 = yc;
        ylc3 = yc;
        
        xlc = [xlc0 xlc1 xlc2 xlc3];
        ylc = [ylc0 ylc1 ylc2 ylc3];
        
%         for j = 2:3
%             
%             con = ylc(j) < ylc(j-1);
%             
%             if con
%                 
%                 ylc(j) = ylc(j-1);
%             end
%             
%             con = ylt(j) < ylt(j-1);
%             
%             if con
%                 
%                 ylt(j) = ylt(j-1);
%             end
%             
%         end
        
        %% 4 degree Bezier curve
        % Trailing thickness curve
        xtt0 = xt;
        % xtt1 = (7*xt - 3*xlt2)/4;
        xtt1 = (7*xt + (9*b8.^2)/(2*rle))/4;
        % xtt2 = 3*xt - 2.5*xlt2;
        xtt2 = 3*xt + (15*b8.^2)/(4*rle);
        xtt3 = b15;
        xtt4 = 1;
        
        ytt0 = yt;
        ytt1 = yt;
        ytt2 = (yt + b8)/2;
        ytt3 = dzte + (1 - b15)*tan(bte);
        ytt4 = dzte;
        
        xtt = [xtt0 xtt1 xtt2 xtt3 xtt4];
        ytt = [ytt0 ytt1 ytt2 ytt3 ytt4];
        
        % Trailing camber curve
        xtc0 = xc;
        % xtc1 = (3*xc - yc*cgle)/2;
        xtc1 = (7*xc - 3*b2)/4;
        % xtc2 = (-8*yc*cgle + 13*xc)/6;
        xtc2 = (b0 - 5*b2 + 6*xc)/2;
        xtc3 = b17;
        xtc4 = 1;
        
        xtc2(xtc2 <= xtc1) = xtc1 + 0.01;
        
        ytc0 = yc;
        ytc1 = yc;
        % ytc2 = 5*yc/6;
        ytc2 = (yc + b1)/2;
        ytc3 = zte - (1 - b17)*tan(ate);
        % ytc3 = zte + (1 - b17)*tan(ate);
        ytc4 = zte;
        
        xtc = [xtc0 xtc1 xtc2 xtc3 xtc4];
        ytc = [ytc0 ytc1 ytc2 ytc3 ytc4];
        
%         for j = 2:4
%             
%             con = ytc(j) > ytc(j-1);
%             
%             if con
%                 
%                 ytc(j) = ytc(j-1);
%             end
%             
%             con = ytt(j) > ytt(j-1);
%             
%             if con
%                 
%                 ytt(j) = ytt(j-1);
%             end
%         end
        
        Pxlt = xlt0*(1 - ult).^3 + 3*xlt1*ult.*(1 - ult).^2 + 3*xlt2*ult.^2.*(1 - ult) + xlt3*ult.^3;
        Pxlc = xlc0*(1 - ulc).^3 + 3*xlc1*ulc.*(1 - ulc).^2 + 3*xlc2*ulc.^2.*(1 - ulc) + xlc3*ulc.^3;
        Pxtt = xtt0*(1 - utt).^4 + 4*xtt1*utt.*(1 - utt).^3 + 6*xtt2*utt.^2.*(1 - utt).^2 + 4*xtt3*utt.^3.*(1 - utt) + xtt4*utt.^4;
        Pxtc = xtc0*(1 - utc).^4 + 4*xtc1*utc.*(1 - utc).^3 + 6*xtc2*utc.^2.*(1 - utc).^2 + 4*xtc3*utc.^3.*(1 - utc) + xtc4*utc.^4;
        
        Pylt = ylt0*(1 - ult).^3 + 3*ylt1*ult.*(1 - ult).^2 + 3*ylt2*ult.^2.*(1 - ult) + ylt3*ult.^3;
        Pylc = ylc0*(1 - ulc).^3 + 3*ylc1*ulc.*(1 - ulc).^2 + 3*ylc2*ulc.^2.*(1 - ulc) + ylc3*ulc.^3;
        Pytt = ytt0*(1 - utt).^4 + 4*ytt1*utt.*(1 - utt).^3 + 6*ytt2*utt.^2.*(1 - utt).^2 + 4*ytt3*utt.^3.*(1 - utt) + ytt4*utt.^4;
        Pytc = ytc0*(1 - utc).^4 + 4*ytc1*utc.*(1 - utc).^3 + 6*ytc2*utc.^2.*(1 - utc).^2 + 4*ytc3*utc.^3.*(1 - utc) + ytc4*utc.^4;
        
        % % One to one, monotonically inc/dec tests
        % xCombineDiff = [diff(Pxlc); diff(Pxlt); diff(Pxtc); diff(Pxtt)];
        % yLeadingDiff = [diff(Pylc); diff(Pylt)];
        % yTrailingDiff = [diff(Pytc); diff(Pytt)];
        % yTrailingDiff = diff(Pytt);
        %
        % horizTest = any(yLeadingDiff < 0) | any(yTrailingDiff > 0);
        % vertiTest = any(xCombineDiff < 0);
        %
        % if horizTest || vertiTest
        %
        %     continue
        % else
        %     true
        % end
        
        Pxc = [Pxlc; Pxtc(2:end)];
        Pyc = [Pylc; Pytc(2:end)];
        Pxt = [Pxlt; Pxtt(2:end)];
        Pyt = [Pylt; Pytt(2:end)];
        
        Pyt = interp1(Pxt, Pyt, Pxc);
        
        %% Camber line feasibility test
        xdiff = diff(Pxc);
        zdiff = diff(Pyc);
        
        xdiff(end+1) = -(Pxc(end-1) - Pxc(end));
        zdiff(end+1) = Pyc(end-1) - Pyc(end);
        
        theta = atan(abs(zdiff)./xdiff);
        
        theta(isnan(theta)) = 0;
        
        con = zdiff < 0;
        theta(con) = -theta(con);
        
        if any(abs(theta) > 50 * pi/180)
            
            flag(i) = true;
            continue
        end
        
        %% Thickness line feasibility test
        xdiff = diff(Pxt);
        zdiff = diff(Pyt);
        
        phi = atan(abs(zdiff)./xdiff);
        
        phi(isnan(phi)) = 0;
        
        con = zdiff < 0;
        phi(con) = -phi(con);
        
        if any(abs(phi(ut + 1:end)) > 50 * pi/180)
            
            flag(i) = true;
            continue
        end
        
        xu = Pxc - Pyt .* sin(theta);% - t .* sin(theta);
        zu = Pyc + Pyt .* cos(theta);% + t .* sin(theta);
        
        xl = Pxc + Pyt .* sin(theta);% - t .* sin(theta);
        zl = Pyc - Pyt .* cos(theta);
        
        % Ensure final values are at x = 1 to avoid NaN in interpolation
        zu(end) = zu(end-1) + (1 - xu(end-1)).*((zu(end) - zu(end-1))./(xu(end) - xu(end-1)));
        zl(end) = zl(end-1) + (1 - xl(end-1)).*((zl(end) - zl(end-1))./(xl(end) - xl(end-1)));
        
        xu(end) = 1;
        xl(end) = 1;
        
        Pzu = interp1(xu, zu, Pxu);
        Pzl = interp1(xl, zl, Pxl);
        
%         figure
%         hold on
%         axis equal
%         plot(Pxlt,Pylt,'r');
%         plot(Pxlc,Pylc,'b');
%         plot(Pxtt,Pytt,'g');
%         plot(Pxtc,Pytc,'y');
%         plot(Pxu,Pzu,'k.');
%         plot(Pxl,Pzl,'k.');
%         plot(xlt,ylt,'rx');
%         plot(xlc,ylc,'bx');
%         plot(xtt,ytt,'gx');
%         plot(xtc,ytc,'yx');
%         xlabel('x/c');
%         ylabel('y/c');
%         title('Aerofoil Section');
        % legend({'Camber','Control Points'},'Location','northeast');
        
        
        if nargin >= 6
            
            % To maintain same format as comparison aerofoil
            foilData{i,j} = [Pxu, Pzu; Pxl(2:end), Pzl(2:end)];
        else
        
            foilData{i,j} = [Pxu, Pzu; flipud([Pxl, Pzl])];
        end
    end
end