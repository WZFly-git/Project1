function  pic_rho_BWC()

dim = 11;
B = linspace(4.6e-6,5.2e-6,dim);    % 4.9087e-06
W = linspace(0.8e-6,1.1e-6,dim);    % 1.0901e-06
C = linspace(3.8e-6,4.5e-6,dim);    % 3.9710e-06
% compute(B, W, C);
% compute_double(B, W, C);
close all

figure('Position', [120,200,600,520]);
data = load('data_next/more_par/cube.mat');
cube_value = data.cube_value;

[X, Y, Z] = meshgrid(B, W, C);
X = X(:);       Y = Y(:);       Z = Z(:);
value = cube_value(:);

sz_min = 10;
sz_max = 100;
dot_size = sz_min + (sz_max - sz_min) * (value - min(value)) / (max(value) - min(value));

scatter3(X, Y, Z, dot_size, value, 'filled');
colormap(autumn);
grid on; box on; hold on;

xlim([0.98*min(B),1.02*max(B)]);      ylim([0.98*min(W),1.02*max(W)]);      zlim([0.98*min(C),1.02*max(C)]);
XT = linspace(4.6e-6,5.2e-6,4);     xticks(XT);
YT = linspace(0.8e-6,1.1e-6,4);     yticks(YT);
ZT = linspace(3.8e-6,4.6e-6,5);     zticks(ZT);
xlabel('\rho_B (nmol/L)');     ylabel('\rho_W (nmol/L)');       zlabel('\rho_C (nmol/L)');
set(gca,'FontSize',12,'FontName','Arial');
print('Figure\rho_BWC\cube_','-dpng','-r600');



BW = load('data_next/double_par/BW_cell.dat');  BW = BW';
WC = load('data_next/double_par/WC_cell.dat');  WC = WC';
CB = load('data_next/double_par/CB_cell.dat');  CB = CB';

sz_min = 15;
sz_max = 150;

figure(2)
[X, Y] = meshgrid(B, W);
X = X(:);       Y = Y(:);   value = BW(:);

dot_size = sz_min + (sz_max - sz_min) * (value - min(value)) / (max(value) - min(value));

scatter(X, Y, dot_size, value, 'filled');
xlim([0.99*min(B),1.01*max(B)]);      ylim([0.99*min(W),1.01*max(W)]);
xlabel('\rho_B (nmol/L)');     ylabel('\rho_W (nmol/L)'); 
colormap(autumn);
grid on; box on; hold on;
colorbar;
set(gca,'FontSize',12,'FontName','Arial');
print('Figure\rho_BWC\double_BW','-dpng','-r600');

figure(3)
[X, Y] = meshgrid(W, C);
X = X(:);       Y = Y(:);   value = WC(:);

dot_size = sz_min + (sz_max - sz_min) * (value - min(value)) / (max(value) - min(value));

scatter(X, Y, dot_size, value, 'filled');
xlim([0.99*min(W),1.01*max(W)]);      ylim([0.99*min(C),1.01*max(C)]);
xlabel('\rho_W (nmol/L)');     ylabel('\rho_C (nmol/L)'); 
colormap(autumn);
grid on; box on; hold on;
colorbar;
set(gca,'FontSize',12,'FontName','Arial');
print('Figure\rho_BWC\double_WC','-dpng','-r600');

figure(4)
[X, Y] = meshgrid(C, B);
X = X(:);       Y = Y(:);   value = CB(:);

dot_size = sz_min + (sz_max - sz_min) * (value - min(value)) / (max(value) - min(value));

scatter(X, Y, dot_size, value, 'filled');
xlim([0.99*min(C),1.01*max(C)]);      ylim([0.99*min(B),1.01*max(B)]);
xlabel('\rho_C (nmol/L)');     ylabel('\rho_B (nmol/L)'); 
colormap(autumn);
grid on; box on; hold on;
colorbar;
set(gca,'FontSize',12,'FontName','Arial');
print('Figure\rho_BWC\double_CB','-dpng','-r600');
end





%% 
function compute(B, W, C)
% define global para and initialization
parameter();
control();

global md CTLA
T = 0:md.h:md.endTime;
y0 = initialization();

% change two-parameter value and calculation
scheme = 3; % anti-CTLA-4 group        
choose_par_durg(scheme);

dim = length(B);
cube_value = zeros(dim,dim,dim);

for i = 1:dim
    CTLA.rho_B = B(i);
    for j = 1:dim
        CTLA.rho_W = W(j);
        for k = 1:dim
            CTLA.rho_C4 = C(k);
            Mat = solve(y0,T);
            cube_value(i,j,k) = Mat(1401,2);
            disp((i-1)*121 + (j-1)*11 + k);
        end
    end
end

save('data_next/more_par/cube.mat', 'cube_value'); 

end


function compute_double(B, W, C)
parameter();
control();

global md CTLA
T = 0:md.h:md.endTime;
y0 = initialization();

% change two-parameter value and calculation
scheme = 3; % anti-CTLA-4 group        
choose_par_durg(scheme);

dim = length(B);
BW_cell = zeros(dim,dim);
WC_cell = zeros(dim,dim);
CB_cell = zeros(dim,dim);

for i = 1:dim
    CTLA.rho_B = B(i);
    for j = 1:dim
        CTLA.rho_W = W(j);
        Mat = solve(y0,T);
        BW_cell(i,j) = Mat(1401,2);
        disp((i-1)*11 + j);
    end
end
parameter();

for i = 1:dim
    CTLA.rho_W = W(i);
    for j = 1:dim
        CTLA.rho_C4 = C(j);
        Mat = solve(y0,T);
        WC_cell(i,j) = Mat(1401,2);
        disp((i-1)*11 + j);
    end
end
parameter();

for i = 1:dim
    CTLA.rho_C4 = C(i);
    for j = 1:dim
        CTLA.rho_B = B(j);
        Mat = solve(y0,T);
        CB_cell(i,j) = Mat(1401,2);
        disp((i-1)*11 + j);
    end
end

writematrix(BW_cell, 'data_next/double_par/BW_cell.dat'); 
writematrix(WC_cell, 'data_next/double_par/WC_cell.dat'); 
writematrix(CB_cell, 'data_next/double_par/CB_cell.dat'); 
end

