function plotmethodsdata(AoA,Mach,methodCN,methodCA,methodCm,impactMatrix,expData)

[rows,~] = size(impactMatrix);

nums = (1:rows)';

labels = cellstr(num2str(nums));
xLim = [AoA(1)-0.5 AoA(end)+0.5]; 

figure
pos = figureposition();
set(gcf, 'Position', pos)
positionVector = [0.35 0.6, 0.35, 0.35];
subplot('Position',positionVector)
hold on
title(['Mach = ' num2str(Mach)]);
grid on
plot(expData(:,1),expData(:,2),'ko')
plot(AoA,methodCN,'k')
for i = 1:rows
    text(AoA,methodCN(i,:),labels{i}, 'VerticalAlignment','bottom',...
        'HorizontalAlignment','right');
end
xlabel(['Angle of Attack (' char(176) ')'])
ylabel('C_N')
xlim(xLim);
hold off

subplot(2,2,3)
% figure
hold on
grid on
plot(expData(:,1),expData(:,3),'ko')
plot(AoA,methodCA,'k')
for i = 1:rows
    text(AoA,methodCA(i,:),labels{i}, 'VerticalAlignment','bottom',...
        'HorizontalAlignment','right');
end
xlabel(['Angle of Attack (' char(176) ')'])
ylabel('C_A')
xlim(xLim);
hold off

subplot(2,2,4)
% subplot(2,2,3)
% figure
hold on
grid on
plot(expData(:,1),expData(:,4),'ko')
plot(AoA,methodCm,'k')
for i = 1:rows
    text(AoA,methodCm(i,:),labels{i}, 'VerticalAlignment','bottom',...
        'HorizontalAlignment','right');
end
xlabel(['Angle of Attack (' char(176) ')'])
ylabel('C_m')
xlim(xLim);
hold off

% figure
% hold on
% grid on
% plot(expData(:,1),expData(:,2),'ko')
% plot(AoA,methodCl,'k')
% for i = 1:rows
%     text(AoA,methodCl(i,:),labels{i}, 'VerticalAlignment','bottom',...
%         'HorizontalAlignment','right');
% end
% xlabel(['Angle of Attack (' char(176) ')'])
% ylabel('C_L')
% hold off
% 
% figure
% hold on
% grid on
% plot(expData(:,1),expData(:,3),'ko')
% plot(AoA,methodCd,'k')
% for i = 1:rows
%     text(AoA,methodCd(i,:),labels{i}, 'VerticalAlignment','bottom',...
%         'HorizontalAlignment','right');
% end
% xlabel(['Angle of Attack (' char(176) ')'])
% ylabel('C_D')
% hold off