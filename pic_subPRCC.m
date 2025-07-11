function pic_subPRCC()

%% prepare and process data
para_mat = load('data/PRCC/para_mat.dat');
target = load('data/PRCC/target.dat');
p = load('data/PRCC/p.dat');

[all_par, fm] = generate();

[p_sort, Index] = sort(p, 'descend');
four_pos = para_mat(:,Index(1:4));
four_neg = para_mat(:,Index(57:60));
par = [four_pos, four_neg];

Y = cell(1, 8);
for s = 1:8
    Y{s} = linear_regression(par(:,s), target);
end


%%
par_nametype = {'$\rho_{C4}$','$T_{N4}$','$\beta_{C}$','$d_{T_c}$',...
               '$\delta_{I_{12} D}$','$\kappa_{T_c}$','$d_{T_\beta}$','$\rho_B$'};

%linecolor = ["#FF4E22","#FFA722","#4EFF22","#227AFF",...
%             "#227AFF","#4EFF22","#FFA722","#FF4E22"];

circlecolor = ['#ea3d33', '#ed803b', '#75eb47', '#55bbeb',...
             "#20a1ed","#4444f5","#a74af6","#ea3395"];
wcirclecolor = circlecolor;
linecolor = circlecolor;



%% plot 
close all
for i = 1:4
    x = par(:,i);
    figure('Position', [i*400-300,480,400,300]);

    scatter(x, target, 40, 'filled', 'MarkerFaceColor', circlecolor{i},...
            'MarkerEdgeColor', wcirclecolor{i}, 'MarkerFaceAlpha', 1) 
    hold on
    box on
    plot(Y{i}(1,:), Y{i}(2,:), 'LineStyle', '-', 'Color', linecolor{i}, 'LineWidth', 1.8);

    xlabel(par_nametype{i}, 'Interpreter', 'latex');
    
    xlim([all_par.(fm{Index(i)})*0.945 all_par.(fm{Index(i)})*1.055]);
    ylim([-0.5 0.4]);

    ax = gca;
    xlims = ax.XLim;

    x_place = xlims(1) + 0.59 * (xlims(2) - xlims(1));
    x_rho = xlims(1) + 0.84 * (xlims(2) - xlims(1));
    y_place = -0.45;
    content1 = sprintf('Sample = %d', length(target));
    content2 = sprintf('Spearman''s    = %.2f', p_sort(i));
    content3 = sprintf('$\\rho $');
    
    text(x_place, y_place+0.08, content1,...
       'FontName', 'Arial', 'FontSize', 10);
    text(x_place, y_place, content2,...
        'FontName', 'Arial', 'FontSize', 10);
    text(x_rho, y_place, content3,...
        'FontName', 'Arial', 'FontSize', 10, 'Interpreter', 'Latex');
    grid on
    set(gca,'FontName','Airal','FontSize',10);
    %print(['Figure/PRCC/sub_',num2str(i)],'-dpng','-r600');
end



for k = 5:8
    x = par(:,k);
    figure('Position', [k*400-1900,100,400,300]);
    
    scatter(x, target, 50, 'filled', 'MarkerFaceColor', circlecolor{k},...
            'MarkerEdgeColor', wcirclecolor{k}, 'MarkerFaceAlpha', 0);
    hold on
    box on
    plot(Y{k}(1,:), Y{k}(2,:), 'LineStyle', ':', 'Color', linecolor{k}, 'LineWidth', 1.8);

    xlabel(par_nametype{k}, 'Interpreter', 'latex');

    xlim([all_par.(fm{Index(52+k)})*0.945 all_par.(fm{Index(52+k)})*1.055]);
    ylim([-0.5 0.4]);

    ax = gca;
    xlims = ax.XLim;

    x_place = xlims(1) + 0.1 * (xlims(2) - xlims(1));
    x_rho = xlims(1) + 0.35 * (xlims(2) - xlims(1));
    y_place = -0.45;

    content1 = sprintf('Sample = %d', length(target));
    content2 = sprintf('Spearman''s    = %.2f', p_sort(k+52));
    content3 = sprintf('$\\rho $');
    
    text(x_place, y_place+0.08, content1,...
       'FontName', 'Arial', 'FontSize', 10);
    text(x_place, y_place, content2,...
        'FontName', 'Arial', 'FontSize', 10);
    text(x_rho, y_place, content3,...
        'FontName', 'Arial', 'FontSize', 10, 'Interpreter', 'Latex');  
    grid on
    set(gca,'FontName','Airal','FontSize',10);
    %print(['Figure/PRCC/sub_',num2str(k)],'-dpng','-r600');
end

end







function [all_par,fm] = generate()

parameter();
control();

global par PDL CTLA Tumor

D1 = 'alpha2_default';      PDL = rmfield(PDL, D1);
D2 = 'K_Q2_default';        CTLA = rmfield(CTLA,D2);
fm1 = fieldnames(par);  fm2 = fieldnames(PDL);  fm3 = fieldnames(CTLA); fm4 = fieldnames(Tumor);
fm = [fm1; fm2; fm3; fm4];
m1 = length(fm1);       m2 = length(fm2);       m3 = length(fm3);   m4 = length(fm4);
M = [m1, m2, m3, m4]; 

all_par = struct();
g_col = 1;
for i = 1:length(M)
    str = choose_str(i);
    for colume = 1:M(i)
        name = fm{g_col};
        all_par.(name) = str.(name); 
        % 10% oscillate
        g_col = g_col + 1;
    end
end

end


function str = choose_str(i)

global par PDL CTLA Tumor

if i == 1
    str = par;
elseif i == 2
    str = PDL;
elseif i == 3
    str = CTLA;
else 
    str = Tumor;
end


end



function Y = linear_regression(x,y)

c = polyfit(x, y, 1);
xa = min(x);    xb = max(x);
n = 50;
x = linspace(xa,xb,n);
y = polyval(c,x);
Y = [x; y];
end

