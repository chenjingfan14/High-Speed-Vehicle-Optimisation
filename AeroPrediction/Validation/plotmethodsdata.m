function plotmethodsdata(AoA,methodCN,methodCA,~,~,methodCm,impactMatrix,~,expData)

[rows,~] = size(impactMatrix);

nums = (1:rows)';

labels = cellstr(num2str(nums));

figure
hold on
grid on
plot(expData(:,1),expData(:,2),'ko')
plot(AoA,methodCN,'k')
for i = 1:rows
    text(AoA,methodCN(i,:),labels{i}, 'VerticalAlignment','bottom',...
        'HorizontalAlignment','right');
end
xlabel('Angle of Attack')
ylabel('C_N')
hold off

figure
hold on
grid on
plot(expData(:,1),expData(:,3),'ko')
plot(AoA,methodCA,'k')
for i = 1:rows
    text(AoA,methodCA(i,:),labels{i}, 'VerticalAlignment','bottom',...
        'HorizontalAlignment','right');
end
xlabel('Angle of Attack')
ylabel('C_A')
hold off

figure
hold on
grid on
plot(expData(:,1),expData(:,4),'ko')
plot(AoA,methodCm,'k')
for i = 1:rows
    text(AoA,methodCm(i,:),labels{i}, 'VerticalAlignment','bottom',...
        'HorizontalAlignment','right');
end
xlabel('Angle of Attack')
ylabel('C_m')
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
% xlabel('Angle of Attack')
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
% xlabel('Angle of Attack')
% ylabel('C_D')
% hold off