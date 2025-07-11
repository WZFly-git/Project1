function pic_Tumor()
% to plot Tumor dynamics under four scheme

vehicle_date  = [100; 270.588; 531.373; 713.725; 1133.333]*8e4;
PDL1_data = [100; 217.647; 421.568; 607.843; 964.705]* 8e4;
CTLA4_data  = [100; 227.451; 372.549; 564.706; 1005.882]* 8e4;
combine_data = [100; 203.922; 278.431; 407.843; 539.216]* 8e4;
real = [vehicle_date, PDL1_data, CTLA4_data, combine_data]';
giveTime = [1,3.5,7,10.5,14];
R = R2();

close all

titletype = {'vehicle','anti-PD-L1','anti-CTLA-4','combine'};
color_vehicle = [0.98,0.106,0.294];
color_anti_PD_L1 = [0.961,0.518,0.039];
color_anti_CTLA_4 = [57,127,199]/255;
color_combine = [0.698,0,0.859];
colortype = {color_vehicle, color_anti_PD_L1, color_anti_CTLA_4, color_combine};

labeltype = {'s','^','o','d'};
name = {'vehicle','anti-PD-L1', 'anti-CTLA-4', 'combine'};

figure('Position', [100,200,1500,300]);
for i=1:4
    filename = ['data/basic/',name{i},'.dat'];
    A = load(filename);
    subplot(1,4,i);
    grid on; box on; hold on;
    plot(A(:,1), A(:,2), 'Color', colortype{i}, 'linewidth', 2);
    scatter(giveTime, real(i,:),100,colortype{i},labeltype{i},'filled');
    % scatter ontrol is sequential
    xlabel('Time (days)');
    ylabel('Tumor Number (cells)')
    xlim([0, 15]); ylim([0, 1.3e8]);
    title(titletype{i});
    text(2, 1.1e8, ['R^2 = ', num2str(round(R(i),2))],'FontSize',12,'Color',colortype{i},'FontWeight','bold');
end

set(gca,'FontName','Airal','FontSize',10);
%print('Figure/Tumor_fit/fit','-dpng','-r600');




figure('Position', [100,200,600,500]);
for i=1:4
    filename = ['data/basic/',name{i},'.dat'];
    A = load(filename);
    grid on; box on; hold on;
    plot(A(:,1), A(:,3), 'Color', colortype{i}, 'linewidth', 2);
    xlabel('Time (days)');
    ylabel('Dendritic Number (cells)')
    xlim([0, 15]); ylim([3e8, 4.5e8]);
end
set(gca,'FontName','Airal','FontSize',18);
print('Figure/Tumor_fit/DCs_','-dpng','-r600');

figure('Position', [100,400,600,500]);
for i=1:4
    filename = ['data/basic/',name{i},'.dat'];
    A = load(filename);
    grid on; box on; hold on;
    plot(A(:,1), A(:,4), 'Color', colortype{i}, 'linewidth', 2);
    xlabel('Time (days)');
    ylabel('Cytotoxic T cell (cells)')
    xlim([0, 15]); ylim([0, 1.2e10]);
end
set(gca,'FontName','Airal','FontSize',18);
print('Figure/Tumor_fit/Tc_','-dpng','-r600');

figure('Position', [100,200,600,500]);
for i=1:4
    filename = ['data/basic/',name{i},'.dat'];
    A = load(filename);
    grid on; box on; hold on;
    plot(A(:,1), A(:,5), 'Color', colortype{i}, 'linewidth', 2);
    xlabel('Time (days)');
    ylabel('Helper T cell (cells)')
    xlim([0, 15]); ylim([0, 3e9]);
end
set(gca,'FontName','Airal','FontSize',18);
print('Figure/Tumor_fit/Th_','-dpng','-r600');

figure('Position', [100,500,600,500]);
for i=1:4
    filename = ['data/basic/',name{i},'.dat'];
    A = load(filename);
    grid on; box on; hold on;
    plot(A(:,1), A(:,6), 'Color', colortype{i}, 'linewidth', 2);
    xlabel('Time (days)');
    ylabel('Regulatory T cell (cells)')
    xlim([0, 15]); ylim([0, 6e9]);
end
set(gca,'FontName','Airal','FontSize',18);
print('Figure/Tumor_fit/Tr_','-dpng','-r600');







% 
% figure('Position', [100,200,1500,300]);
% for i=1:4
%     filename = ['data/basic/',name{i},'.dat'];
%     A = load(filename);
%     subplot(1,4,i);
%     grid on; box on; hold on;
%     plot(A(:,1), A(:,3), 'Color', colortype{i}, 'linewidth', 2);
%     xlabel('Time (days)');
%     ylabel('Dendritic Number (cells)')
%     xlim([0, 15]); ylim([3e8, 4.5e8]);
% end
% set(gca,'FontName','Airal','FontSize',10);
% print('Figure/Tumor_fit/DCs','-dpng','-r600');
% 
% figure('Position', [100,400,1500,300]);
% for i=1:4
%     filename = ['data/basic/',name{i},'.dat'];
%     A = load(filename);
%     subplot(1,4,i);
%     grid on; box on; hold on;
%     plot(A(:,1), A(:,4), 'Color', colortype{i}, 'linewidth', 2);
%     xlabel('Time (days)');
%     ylabel('Cytotoxic T cell (cells)')
%     xlim([0, 15]); ylim([0, 1.2e10]);
% end
% set(gca,'FontName','Airal','FontSize',10);
% print('Figure/Tumor_fit/Tc','-dpng','-r600');
% 
% figure('Position', [100,200,1500,300]);
% for i=1:4
%     filename = ['data/basic/',name{i},'.dat'];
%     A = load(filename);
%     subplot(1,4,i);
%     grid on; box on; hold on;
%     plot(A(:,1), A(:,5), 'Color', colortype{i}, 'linewidth', 2);
%     xlabel('Time (days)');
%     ylabel('Helper T cell (cells)')
%     xlim([0, 15]); ylim([0, 3e9]);
% end
% set(gca,'FontName','Airal','FontSize',10);
% print('Figure/Tumor_fit/Th','-dpng','-r600');
% 
% figure('Position', [100,500,1500,300]);
% for i=1:4
%     filename = ['data/basic/',name{i},'.dat'];
%     A = load(filename);
%     subplot(1,4,i);
%     grid on; box on; hold on;
%     plot(A(:,1), A(:,6), 'Color', colortype{i}, 'linewidth', 2);
%     xlabel('Time (days)');
%     ylabel('Regulatory T cell (cells)')
%     xlim([0, 15]); ylim([0, 6e9]);
% end
% set(gca,'FontName','Airal','FontSize',10);
% print('Figure/Tumor_fit/Tr','-dpng','-r600');

end