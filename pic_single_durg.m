function pic_single_durg()

data = load('data_next/Durg/single_Q1.mat');
A = data.A;
T = 0:0.01:14;
% plot figure
dim = 11;
Q1 = linspace(0,2,dim); 

linestype={'k-','b','m-','y-','r-','k'};
cellstype = {'Tumor cell','Cytotoxic T cell','Helper T cell','Regulatory T cell'};
close all

% change rho_P 
figure('Position', [120,700,1500,300]);
for k = 1:6
    color = linestype{k};
    for i=1:4
        subplot(1,4,i);
        grid on; box on; hold on;
        plot(T,A{i}(:,k*2-1),color,'linewidth',2);
        xlabel('Time(days)');
        title(cellstype{i});
        hold on
    end
end

subplot(1,4,1);
leg = cell(1,6);
for k = 1:6
    leg{k} = sprintf('\\rho_P = %f', Q1(k*2-1));
end
legend(leg, 'Location', 'northwest','FontSize',8);
legend('boxoff');

end




function com()
parameter();
control();

global md Dose
T = 0:md.h:md.endTime;
y0 = initialization();

scheme = 2; % combine group
choose_par_durg(scheme);

dim = 11; % matrix dimension
A = cell(1,5);
for i = 1:5
    A{i} = zeros(1401,dim);
end
Q2 = linspace(0,2,dim); 

%% calculation and save data
for s = 1:dim
    Dose.Q2 = Q2(s);
    Mat = solve(y0, T);
    for i = 1:5
        A{i}(:,s) = Mat(1:1401, i+1);    
    end
end

save('data_next/Durg/single_Q2.mat','A');


end




























