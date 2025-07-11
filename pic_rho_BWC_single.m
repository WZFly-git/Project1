function pic_rho_BWC_single()

dim = 11;
B = linspace(4.6e-6,5.2e-6,dim);    % 4.9087e-06
W = linspace(0.8e-6,1.1e-6,dim);    % 1.0901e-06
C = linspace(3.8e-6,4.5e-6,dim);    % 3.9710e-06
% compute(B,W,C);
close all

T = 0:0.01:14;
colorstype = {'#ff0c00', '#ff8c00', '#ffd700', '#00ff00', '#4169e1', '#9400d3'};
cellstype = {'Tumor cell','Cytotoxic T cell','Helper T cell','Regulatory T cell'};

figure('Position', [120,50,1500,200]);

data = load('data_next/single_par/B_cell.mat');
Number = data.Number_B;   Number(2) = [];
for k = 1:6
    color = colorstype{k};
    for i=1:4
        subplot(1,4,i);
        grid on; box on; hold on;
        plot(T,Number{i}(:,k*2-1),'Color', color, 'linewidth', 1.2);
        xlabel('Time (days)');
        ylabel(cellstype{i})
        hold on
    end
end

subplot(1,4,1);
leg = cell(1,6);
for k = 1:6
    leg{k} = sprintf('\\rho_B = %.2e', B(k*2-1));
end
legend(leg, 'Location', 'northwest','FontSize',6);
legend('boxoff');
print('Figure\rho_BWC\singleB','-dpng','-r600');



figure('Position', [120,600,1500,200]);

data = load('data_next/single_par/W_cell.mat');
Number = data.Number_W;   Number(2) = [];
for k = 1:6
    color = colorstype{k};
    for i=1:4
        subplot(1,4,i);
        grid on; box on; hold on;
        plot(T,Number{i}(:,k*2-1),'Color', color, 'linewidth', 1.2);
        xlabel('Time (days)');
        ylabel(cellstype{i})
        hold on
    end
end

subplot(1,4,1);
leg = cell(1,6);
for k = 1:6
    leg{k} = sprintf('\\rho_W = %.2e', W(k*2-1));
end
legend(leg, 'Location', 'northwest','FontSize',6);
legend('boxoff');
print('Figure\rho_BWC\singleW','-dpng','-r600');



figure('Position', [120,300,1500,200]);

data = load('data_next/single_par/C_cell.mat');
Number = data.Number_C;   Number(2) = [];
for k = 1:6
    color = colorstype{k};
    for i=1:4
        subplot(1,4,i);
        grid on; box on; hold on;
        plot(T,Number{i}(:,k*2-1),'Color', color, 'linewidth', 1.2);
        xlabel('Time (days)');
        ylabel(cellstype{i})
        hold on
    end
end

subplot(1,4,1);
leg = cell(1,6);
for k = 1:6
    leg{k} = sprintf('\\rho_{C_4} = %.2e', C(k*2-1));
end
legend(leg, 'Location', 'northwest','FontSize',6);
legend('boxoff');
print('Figure\rho_BWC\singleC','-dpng','-r600');




end


function compute(B,W,C)
% define global para and initialization
parameter();
control();

global md CTLA
T = 0:md.h:md.endTime;
y0 = initialization();

% change two-parameter value and calculation
scheme = 3; % anti-PD-L1 group        
choose_par_durg(scheme);

Number_B = cell(1,5);     Number_W = cell(1,5);     Number_C = cell(1,5);
dim = length(B);
for k = 1:5
    Number_B{k} = zeros(1401, dim);
    Number_W{k} = zeros(1401, dim);
    Number_C{k} = zeros(1401, dim);
end



for s = 1:dim
    CTLA.rho_B = B(s);
    Mat = solve(y0,T);
    for k = 1:5
        Number_B{k}(:, s) = Mat(1:1401, k+1);
    end
end

parameter();

for t = 1:dim
    CTLA.rho_W = W(t);
    Mat = solve(y0,T);
    for k = 1:5
        Number_W{k}(:, t) = Mat(1:1401, k+1);
    end
end

parameter();

for r = 1:dim
    CTLA.rho_C4 = C(r);
    Mat = solve(y0,T);
    for k = 1:5
        Number_C{k}(:, r) = Mat(1:1401, k+1);
    end
end

save('data_next/single_par/B_cell.mat', 'Number_B'); 
save('data_next/single_par/W_cell.mat', 'Number_W'); 
save('data_next/single_par/C_cell.mat', 'Number_C'); 

end



