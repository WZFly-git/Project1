function pic_double_durg()
% Curves under different parameters --- rho_P and rho_L
% combine scheme

%double_dose();
%double_dose_bestplan();


dim = 11;
close all

%% 
figure('Position', [50,200,800,600]);
M = load('data_next/Durg/double_durg.dat'); 
M = M';
x = linspace(0, 10, dim);
y = linspace(0, 10, dim);
[X, Y] = meshgrid(x, y);
pcolor(X, Y, M);      % auto turn matrix
shading interp
hold on; box on; grid on;
contour(X, Y, M, 10, 'LineWidth', 1, 'Color', 'k', 'LineStyle', ':');  
plot(9.9, 9.9, 'p', 'MarkerSize', 6, 'LineWidth', 2,'MarkerFaceColor', 'w','Color', 'w');


xlabel('anti-PD-L1 (nmol)');        ylabel('anti-CTLA-4 (nmol)');      
xticks(0:2:10);        yticks(0:2:10);
%title('The combined treatment effect under different Dose');
colormap(jet);
colorbar
set(gca,'FontSize',14,'FontName','Arial');
% print('Figure/Durg/dose_new','-dpng','-r600');


% %%
figure('Position', [100,300,800,600]);
M_best = load('data_next/Durg/double_durg_bestplan.dat');
M_best = M_best';
x = linspace(0, 10, dim); 
y = linspace(0, 10, dim); 
[X, Y] = meshgrid(x, y);
pcolor(X, Y, M_best);      % auto turn matrix
shading interp
hold on; box on; grid on;
contour(X, Y, M_best, 5, 'LineWidth', 1, 'Color', 'k', 'LineStyle', ':');  
contour(X, Y, M_best, [M(end,end), M(end,end)], 'LineWidth', 2, 'Color', 'w', 'LineStyle', ':');  

xlabel('anti-PD-L1 (nmol)');    ylabel('anti-CTLA-4 (nmol)');
xticks(0:2.5:10);        yticks(0:2.5:10);
% title('bestplan');
colormap(jet);
colorbar
set(gca,'FontSize',14,'FontName','Arial');
%print('Figure/Durg/dose_bestplan','-dpng','-r600');








%% 
figure('Position', [900,200,800,600]);
x_fit = linspace(0, 10, 11); 
y_fit = linspace(0, 10, 11); 
Z = fit_data(x_fit,y_fit);
[X_, Y_] = meshgrid(x_fit, y_fit);
surf(X_,Y_,Z);
shading interp
hold on; grid on; box on; 

xlabel('anti-PD-L1 (nmol)');        ylabel('anti-CTLA-4 (nmol)');       zlabel('e(\gamma_A, \gamma_B)');
xticks(0:2:10);        yticks(0:2:10);        zlim([0,1]);
colormap(jet);
colorbar
set(gca,'FontSize',14,'FontName','Arial');
print('Figure/Durg/fit_new','-dpng','-r600');



%%
figure(4);
M = M(:);   Z = Z(:);

scatter(M, Z, 18, 'filled', 'MarkerFaceColor', [0.8,0.8,0.8],...
            'MarkerEdgeColor', [0.8,0.8,0.8], 'MarkerFaceAlpha', 1);
hold on;  box on;
plot(0:0.01:1, 0:0.01:1, 'LineStyle', ':', 'Color', [0.85,0.85,0.85], 'LineWidth', 1);
xlim([0 1]);    ylim([0 1]);
xticks(0:0.2:1);    yticks(0:0.2:1);
xlabel('E(\gamma_A, \gamma_B)');        ylabel('e(\gamma_A, \gamma_B)');
set(gca,'FontSize',10,'FontName','Arial');
print('Figure/Durg/deviation_new','-dpng','-r600');





%%

end










function double_dose()

parameter();
control();

global md 
T = 0:md.h:md.endTime_dose;
y0 = initialization();

scheme = 4; % combine group
choose_par_durg(scheme);

dim = 11; % matrix dimension

Q1 = linspace(0,2,dim); 
Q2 = linspace(0,2,dim); 
A = zeros(dim);

%% calculation and save data
for s = 1:dim
    for t = 1:dim
        md.G1{2} = [1, 1, 1, 1, 1]* Q1(s);
        md.G2{2} = [1, 1, 1, 1, 1]* Q2(t);
        Mat = solve(y0,T);
        A(s,t) = Mat(2101,2);
        disp(s+t);
    end
end

star = A(1,1);
B = zeros(dim);
for i = 1:dim
    for j = 1:dim
        B(i,j) = ( star - A(i,j) ) / star;
    end
end
writematrix(B,'data_next/Durg/double_durg.dat');

end



function Z = fit_data(x_fit, y_fit)
B = load('data_next/Durg/double_durg.dat');
B = B';
dim = length(B);
x = 0:1:10;    y = 0:1:10;
x = x';         y = y';
A = zeros(dim*dim,5);
b = zeros(dim*dim,1);

for i = 1:dim
    A(i*11-10:i*11,:) = [ones(dim,1)*x(i), y, ones(dim,1)*(x(i)^2), x(i).*y, y.^2];
    b(i*11-10:i*11) = B(i,:)';
end

coff = A\b;

%coff = round(coff, 2);
for i = 1:length(coff)
    ab{i} = sprintf('%.2e', coff(i));
    coff(i) = str2num(ab{i});
end


m = length(x_fit);      n = length(y_fit);
Z = zeros(m,n);
for i = 1:m
    for j = 1:n
        Z(i,j) = [x_fit(i), y_fit(j), x_fit(i)^2, x_fit(i)*y_fit(j), y_fit(j)^2] * coff;
    end
end




writematrix(coff,'data_next/Durg/fit_coff.dat');

end


function double_dose_bestplan()

parameter();
control();

global md 
T = 0:md.h:md.endTime_dose;
y0 = initialization();

scheme = 4; % combine group
choose_par_durg(scheme);

dim = 11; % matrix dimension

Q1 = linspace(0,2,dim); 
Q2 = linspace(0,2,dim); 
A = zeros(dim);

%% calculation and save data

for s = 1:dim
    for t = 1:dim
        md.G1 = {0, Q1(s)*5};
        md.G2 = {0, Q2(t)*5};
        Mat = solve(y0,T);
        A(s,t) = Mat(2101,2);
        disp(s+t);
    end
end

star = A(1,1);
B = zeros(dim);
for i = 1:dim
    for j = 1:dim
        B(i,j) = ( star - A(i,j) ) / star;
    end
end
writematrix(B,'data_next/Durg/double_durg_bestplan.dat');

end