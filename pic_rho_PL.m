function pic_rho_PL()

%% cancel annotation when need calculation 

% double_rho_PL(); 
% single_rho_PL();

%to_plot_double();
to_plot_single();

end 



%% Draw the combined effect of the two-parameter changes
function to_plot_double()

data = load('data_next/double_par/PL_cell.mat');
Relative_change = data.Relative_change;

% plot figure

P = logspace(-7,-6,11);
L = linspace(1e-6,2e-6,11);

Para = load('input/key_para.dat');
rho_P = Para(5);        rho_L = Para(6);

[X, Y] = meshgrid(L, P);
cellstype = {'Tumor cell','Dendritic cell','Cytotoxic T cell','Helper T cell','Regulatory T cell'};
filename = {'Tumor', 'DC', 'Tc', 'Th', 'Tr'};
close all

for k = 1:5
    figure('Position', [k*250,300,600,450]);
    pcolor(X, Y, -Relative_change{k});      % auto turn matrix
    set(gca, 'YScale', 'log');
    shading interp
    hold on
    contour(X, Y, -Relative_change{k}, [0 0], '--w','LineWidth', 2, 'Color', 'w');
    plot(rho_L, rho_P, 'p', 'MarkerSize', 10, 'LineWidth', 2,'MarkerFaceColor', 'w' ,'Color', 'w');
    xlabel('\rho_L (nmol/L)');      ylabel('\rho_P (nmol/L)');
    title(cellstype{k});
    colormap(jet);
    colorbar
    set(gca,'FontSize',18,'FontName','Arial');
    %print(['Figure\rho_PL\double__',filename{k}],'-dpng','-r600'); % -r is resolution.
end


end


%% Draw the kinetic curves of each component when a single parameter changes
function to_plot_single()

data1 = load('data_next/single_par/P_cell.mat');
Number_P = data1.Number_P;
Number_P(2) = []; % delete DCs colume 
data2 = load('data_next/single_par/L_cell.mat');
Number_L = data2.Number_L;
Number_L(2) = [];
T = 0:0.01:14;
% plot figure

P = logspace(-7,-6,11);
L = linspace(1e-6,2e-6,11);

% colorstype = {[1, 0.388, 0.278], [1, 0.549, 0], [1, 0.843, 0], [0, 1, 0],...
%               [0.255, 0.412, 0.882], [0.58, 0, 0.827]};
colorstype={'#ff0c00', '#ff8c00', '#ffd700', '#00ff00', '#4169e1', '#9400d3'};

cellstype = {'Tumor cell','Cytotoxic T cell','Helper T cell','Regulatory T cell'};
close all

% change rho_P 
figure('Position', [120,200,1500,200]);
for k = 1:6
    color = colorstype{k};
    for i=1:4
        subplot(1,4,i);
        grid on; box on; hold on;
        plot(T, Number_P{i}(:,k*2-1), 'Color', color, 'linewidth', 1.2);
        xlabel('Time (days)');
        ylabel(cellstype{i});
        hold on
    end
end

subplot(1,4,1);
leg_P = cell(1,6);
for k = 1:6
    leg_P{k} = sprintf('\\rho_P = %.2e', P(k*2-1));
end
legend(leg_P, 'Location', 'northwest','FontSize',6);
legend('boxoff');
print('Figure\rho_PL\single_P','-dpng','-r600'); % -r is resolution.


% change rho_L
figure('Position', [120,500,1500,200]);
for k = 1:6
    color = colorstype{k};
    for i=1:4
        subplot(1,4,i);
        grid on; box on; hold on;
        plot(T, Number_L{i}(:,k*2-1), 'Color', color, 'linewidth', 1.2);
        xlabel('Time (days)');
        ylabel(cellstype{i});
        hold on
    end
end

subplot(1,4,1);
leg_L = cell(1,6);
for k = 1:6
    leg_L{k} = sprintf('\\rho_L = %.2e', L(k*2-1));
end
legend(leg_L, 'Location', 'northwest','FontSize',6);
legend('boxoff');
print('Figure\rho_PL\single_L','-dpng','-r600'); % -r is resolution.

end





%% two-parametric analysis -- rho_P and rho_L
function double_rho_PL()
% define global para and initialization
parameter();
control();

global md PDL
T = 0:md.h:md.endTime;
y0 = initialization();

% load basic data
basic = load('data/basic/anti-PD-L1.dat');
star = basic(1401,2:6);

% change two-parameter value and calculation
scheme = 2; % anti-PD-L1 group        
choose_par_durg(scheme);

dim = 11; % matrix dimension
P = logspace(-7,-6,dim);
L = linspace(1e-6,2e-6,dim);
Relative_change = cell(1,5);
for k = 1:5
    Relative_change{k} = zeros(dim);
end

for s = 1:dim
    for t = 1:dim
        PDL.rho_P = P(s);
        PDL.rho_L = L(t);
        Mat = solve(y0,T);
        
        for k = 1:5
            Relative_change{k}(s,t) = (star(k) - Mat(1401,k+1)) / star(k);
        end
    end
end

save('data_next/double_par/PL_cell.mat', 'Relative_change'); 
% if hope save celldata, second control value need ''
% dat file need trird control value: '-binary'   mean: 二进制

end




%% single-parametric analysis -- rho_P and rho_L
function single_rho_PL()
% define global para and initialization
parameter();
control();

global md PDL
T = 0:md.h:md.endTime;
y0 = initialization();

% change two-parameter value and calculation
scheme = 2; % anti-PD-L1 group        
choose_par_durg(scheme);

dim = 11;
P = logspace(-7,-6,dim);
L = linspace(1e-6,2e-6,dim);
Number_P = cell(1,5);
Number_L = cell(1,5);

for k = 1:5
    Number_P{k} = zeros(1401, dim);
    Number_L{k} = zeros(1401, dim);
end

for s = 1:dim
    PDL.rho_P = P(s);
    Mat = solve(y0,T);
    for k = 1:5
        Number_P{k}(:, s) = Mat(1:1401, k+1);
    end
end

parameter();

for t = 1:dim
    PDL.rho_L = L(t);
    Mat = solve(y0,T);
    for k = 1:5
        Number_L{k}(:, t) = Mat(1:1401, k+1);
    end
end

save('data_next/single_par/P_cell.mat', 'Number_P'); 
save('data_next/single_par/L_cell.mat', 'Number_L'); 
% if hope save celldata, second control value need ''
% dat file need trird control value: '-binary'   mean: 二进制


end