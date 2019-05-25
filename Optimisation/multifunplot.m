function multifunplot(it,nonDomFitDisp,options,nFun,limits,fcount)

if options.Baseline
    
    baselineCost = options.Base.Results.Cost;
else
    baselineCost = zeros(1,nFun);
end

labels = options.CostFunctionLabels;

% Have max limits slightly above that of PF so that points are not
% lying on the edges
maxf = max([nonDomFitDisp; baselineCost] ,[],1)*1.2;
minf = min([nonDomFitDisp; baselineCost] ,[],1);
figure(fcount)
clf
title(['Pareto Front (Iteration: ' num2str(it) ')'])
set(gcf, 'Position', [0, 0, 1920, 1200])
hold on
grid on
xlabel(['$\bar{' labels{1} '}$'],'Interpreter','Latex')
ylabel(['$\bar{' labels{2} '}$'],'Interpreter','Latex')

% Max amount of cost functions that can be plotted is 3, so if
% there are > 3, only plot first 3
if nFun > 3
    maxf = maxf(1:3);
    minf = minf(1:3);
end

for i = 1:nFun
    if inv(i) == 0
        limits(i,:) = [min(0,minf(i)), maxf(i)];
    else
        limits(i,:) = [minf(i), maxf(i)];
    end
end

xlim(limits(1,:))
ylim(limits(2,:))

if nFun == 2
    
    plot(baselineCost(:,1), baselineCost(:,2),'ko');
    % plot(parFitDisp(:,1), parFitDisp(:,2),'k*');
    plot(nonDomFitDisp(:,1), nonDomFitDisp(:,2),'kx');
else
    plot3(baselineCost(:,1), baselineCost(:,2), baselineCost(:,3),'ko');
    % plot3(parFitDisp(:,1), parFitDisp(:,2), parFitDisp(:,3),'k*');
    plot3(nonDomFitDisp(:,1), nonDomFitDisp(:,2), nonDomFitDisp(:,3),'kx');
    zlim(limits(3,:))
    zlabel(['$\bar{' labels{3} '}$'],'Interpreter','Latex')
end

legend('Baseline', 'Pareto Front Design')
hold off