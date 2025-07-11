function pic_Bliss()

Matrix = load('data_next/Durg/Bliss.dat');
dim = 10;
Q1 = linspace(0.2,2,dim); 
Q2 = linspace(0.2,2,dim); 

[X, Y] = meshgrid(Q1, Q2);
close all


figure('Position', [250,300,600,450]);
pcolor(X, Y, Matrix);      % auto turn matrix
shading interp
hold on
xlabel('anti-PD-L1 (nmol)');      ylabel('anti-CTLA-4 (nmol)');

colormap(jet);
colorbar
set(gca,'FontSize',18,'FontName','Arial');
%print('Figure/Durg/Bliss','-dpng','-r600'); 






end










function Bliss()

parameter();
control();

global md 
T = 0:md.h:md.endTime_dose;
y0 = initialization();

dim = 10; % matrix dimension

Q1 = linspace(0.2,2,dim); 
Q2 = linspace(0.2,2,dim); 
A = cell(1,3);
A{1} = zeros(1,dim);
A{2} = zeros(1,dim);
A{3} = zeros(dim,dim);


%% calculation and save data
% star
scheme = 1;
choose_par_durg(scheme);
Mat = solve(y0,T);
star = Mat(2101,2);

% anti-PD-L1
scheme = 2;
choose_par_durg(scheme);

for s = 1:dim
    md.G1{2} = [1, 1, 1, 1, 1]* Q1(s);
    Mat = solve(y0,T);
    A{1}(s) = Mat(2101,2);
    disp(s)
end

% anti-CTLA-4
parameter();
control();

scheme = 3;
choose_par_durg(scheme);

for t = 1:dim
    md.G2{2} = [1, 1, 1, 1, 1]* Q2(t);
    Mat = solve(y0,T);
    A{2}(t) = Mat(2101,2);
    disp(t)
end

% combine
parameter();
control();

scheme = 4;
choose_par_durg(scheme);

for s = 1:dim
    for t = 1:dim
        md.G1{2} = [1, 1, 1, 1, 1]* Q1(s);
        md.G2{2} = [1, 1, 1, 1, 1]* Q2(t);
        Mat = solve(y0,T);
        A{3}(s,t) = Mat(2101,2);
        disp(s+t);
    end
end


G1 = zeros(1,dim);
G2 = zeros(1,dim);
G3 = zeros(dim,dim);

for i = 1:dim
    G1(i) = ( star - A{1}(i) ) / star;
    G2(i) = ( star - A{2}(i) ) / star;    
    for j = 1:dim
        G3(i,j) = ( star - A{3}(i,j) ) / star;
    end
end

CI = zeros(dim);
for i = 1:dim
    for j = 1:dim
        frac = G1(i) + G2(j) - G1(i)* G2(j);
        CI(i,j) = frac / G3(i,j);
    end
end

writematrix(CI, 'data_next/Durg/Bliss.dat');


end