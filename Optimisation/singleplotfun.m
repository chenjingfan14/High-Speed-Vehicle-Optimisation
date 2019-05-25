function singleplotfun(it,GlobalBestFitDisp,options,fcount)

labels = options.CostFunctionLabels;

figure(fcount)
clf
title(['Cost Function Time History (Iteration: ' num2str(it) ')'])
set(gcf, 'Position', [0, 0, 1920, 1200])
hold on
grid on
xlabel('Iterations')
ylabel(['$\bar{' labels{1} '}$'],'Interpreter','Latex')

plot(0:it, GlobalBestFitDisp(1:it-1),'k');

if options.Baseline
    
    plot(it, options.Base.Results.Cost,'kx');
    legend('Optimal Design', 'Baseline')
else
    legend('Optimal Design')
end
hold off