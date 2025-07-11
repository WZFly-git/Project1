function pic_PRCC()
%%
% to drow global sensitivity analysis of all parameters
% PRCC based on the number of "tumor cells on the last day"
% use combine group

parameter();
control();

global par PDL CTLA Tumor

D1 = 'alpha2_default';      PDL = rmfield(PDL, D1);
D2 = 'K_Q2_default';        CTLA = rmfield(CTLA,D2);
fm1 = fieldnames(par);  fm2 = fieldnames(PDL);  fm3 = fieldnames(CTLA); fm4 = fieldnames(Tumor);
fm = [fm1; fm2; fm3; fm4];

m1 = length(fm1);       m2 = length(fm2);       m3 = length(fm3);   m4 = length(fm4);
M = [m1, m2, m3, m4];       

%p = PRCC(100, fm, M);
p = load('data/PRCC/p.dat');
[p_sort, Index] = sort(p, 'descend');

namestype = {'$K_C$','$K_D$','$K_{I_2}$','$K_{I_{12}}$','$K_{T_\beta}$','$K_{T_c I_{10}}$','$K_{T_c T_\beta}$',...
    '$K_{T_h I_{10}}$','$K_{T_h T_\beta}$','$K_{T_r I_\gamma}$','$D_0$','$T_{N8}$','$T_{N4}$','$\beta_{T_c I_2}$',...
    '$\beta_{T_h I_2}$','$\lambda_{DC}$','$d_D$','$d_{T_c}$','$d_{T_h}$','$d_{T_r}$','$\kappa_{T_c}$','$\kappa_{T_h}$',...
    '$\kappa_{T_r}$','$\delta_{I_2 T_h}$','$d_{I_2}$','$\delta_{I_{10} T_r}$','$\delta_{I_{10} C}$','$d_{I_10}$',...
    '$\delta_{I_{12} D}$','$d_{I_{12}}$','$\delta_{T_\beta T_r}$','$\delta_{T_\beta C}$','$d_{T_\beta}$',...
    '$\delta_{I_\gamma T_h}$','$\delta_{I_\gamma T_c}$','$d_{I_\gamma}$','$n$','$\varepsilon_C$','$\varepsilon_r$',...
    '$\rho_P$', '$\rho_L$','$K_PL$','$\alpha_1$','$\alpha_2$','$F_1$','$V_1$','$k_{e1}$','$\rho_B$','$\rho_{C4}$','$\rho_W$','$K_{WB}$',... 
    '$K_{Q2}$','$K_{C4}$','$F_2$','$V_2$','$k_{e2}$','$\beta_C$','$C_M$','$\eta_{T_c}$','$\eta_{T_h}$'};
namestype = namestype(Index);

bars = bar(p_sort);
colormap(jet);
bars.CData = p_sort;
bars.FaceColor = 'flat';

caxis([min(p_sort)-0.1, max(p_sort)+0.1]); % very nice function

% colorbar;


grid on;
xticks(1:60);
xticklabels(namestype);
set(gca,'TickLabelInterpreter','latex');
title('Global sensitivity analysis of all parameters');
set(gcf,'unit','centimeters','position',[1 2 40 10]);

%print('Figure/PRCC/main','-dpng','-r600');
end











%% calculation all parameter PRCC value
function p_value = PRCC(n, fm, M)
% prepare data
global md 
y0 = initialization();
T = 0:md.h:md.endTime;

% parameter sampling
m = sum(M);
sample_mat = lhsdesign(n, m);
%sample_mat = rand(n, m);
para_mat = def_par_matrix(sample_mat, fm, M);
target = calculation(n, para_mat, fm, M, y0, T);

% to rank matrix
for colume = 1:m
    para_mat(:, colume) = tiedrank(para_mat(:, colume));
end
outcome = tiedrank(target);

% get PRCC value 
p_value = zeros(1,m);

% Partial correlation coefficient
% for k = 1:m
%     col = setdiff(1:m, k);
%     p_value(k) = partialcorr(outcome, para_mat(:,k), para_mat(:,col));
% end

% Spearman's rank correlation coefficient
for k = 1:m
    p_value(k) = 1 - 6*sum((para_mat(:,k) - outcome).^2) / (n*(n^2-1));
end


writematrix(p_value,'data/PRCC/p.dat');
end



%% parameter sampling by LHS
function para_mat = def_par_matrix(sample_mat, fm, M)

para_mat = zeros(size(sample_mat));

g_col = 1;
for i = 1:length(M)
    str = choose_str(i);

    for colume = 1:M(i)
        name = fm{g_col};
        para_mat(:, g_col) = str.(name)* (0.95 + 0.1*sample_mat(:, g_col)); 
        % 5% oscillate
        g_col = g_col + 1;
    end
end
writematrix(para_mat,'data/PRCC/para_mat.dat');

end









%% Calculate the indicator value
function target = calculation(n, para_mat, fm, M, y0, T)

% basic = load('data/basic/combine.dat');
% bottom = trapz(basic(1:1401,1),basic(1:1401,2));

global par PDL CTLA Tumor
parameter();
target = zeros(n,1);
for s = 1:n
    g_col = 1;
    for i = 1:length(M)
        for colume = 1:M(i)
            name = fm{g_col};
            switch i
                case 1
                    par.(name) = para_mat(s, g_col);
                case 2
                    PDL.(name) = para_mat(s, g_col);
                case 3
                    CTLA.(name) = para_mat(s, g_col);
                case 4
                    Tumor.(name) = para_mat(s, g_col);
            end
            g_col = g_col + 1;
        end
    end
    choose_par_durg(4);
    Mat = solve(y0, T);
    A = polyfit(Mat(1:1401,1), log(Mat(1:1401,2)), 1);
    target(s) = A(1);
    disp(s);
end

writematrix(target,'data/PRCC/target.dat');

end



%% def parameter struct copy
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