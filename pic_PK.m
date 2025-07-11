function pic_PK()

parameter();
control();

global md 
T = 0:md.h:md.endTime_dose;
n = length(T);
Q1 = zeros(1, n);   Q2 = zeros(1, n);


Plan1 = md.Q1_Plan;
Plan2 = md.Q2_Plan;
close all

%% PK cover
figure('Position', [0,0,2000,700]);
colors = {'#ea3d33', '#ed803b', '#75eb47', '#55bbeb', "#ea3395", "#4444f5"};
for r = 1:length(Plan1)
    md.G1 = Plan1{r};
    md.G2 = Plan2{r};
    subplot(2,3,r);
    for l = 1:n
        Q1(l) = Durg1(T(l));
        Q2(l) = Durg2(T(l));
    end
    plot(T, Q1, 'Color', colors{r}, 'LineWidth', 1);
    grid on; box on; hold on;
    plot(T, Q2, '--', 'Color', colors{r}, 'LineWidth', 1);
    %ylim([0 3])
    xlabel('Time (day)');
    ylabel('Durg concentration (nmol/L)');
    set(gca,'FontSize',12,'FontName','Arial');
end
%print('Figure/durg_plan/PK', '-dpng', '-r600');

%% therapeutic effect
figure('Position', [1100,100,600,300]);
target_1 = load('data/Durg_Plan/Q1_target.dat');
target_2 = load('data/Durg_Plan/Q2_target.dat');
data = [target_1, target_2];
b = bar(data);
grid on; box on; hold on;
colormap(parula);
colorbar;
xticks(1:6);
xticklabels({'plan1', 'plan2', 'plan3', 'plan4', 'plan5', 'plan6',});
set(b(1), 'CData', target_1, 'FaceColor', 'flat');
set(b(2), 'CData', target_2, 'FaceColor', 'flat');
caxis([min(target_1)-0.1, max(target_1)+0.1]); % magical code
set(gca,'FontSize',12,'FontName','Arial');
%print('Figure/durg_plan/target', '-dpng', '-r600');



end









function Q = Durg1(t)
% anti-PD-L1
global md PDL 
G = md.G1;
T_Q1 = G{1};
Dose_Q1 = G{2};
Q = 0;  
for i = 1:length(T_Q1)
    if t >= T_Q1(i)
        Q = Q + PDL.F_1* Dose_Q1(i)* exp(-PDL.k_e1* (t - T_Q1(i)) )/PDL.V_1;
    end
end

end

function Q = Durg2(t)
% anti-CTLA-4
global md CTLA
G = md.G2;
T_Q2 = G{1};
Dose_Q2 = G{2};

Q = 0;  
for i = 1:length(T_Q2)
    if t >= T_Q2(i)
        Q = Q + CTLA.F_2* Dose_Q2(i)* exp(-CTLA.k_e2* (t - T_Q2(i)) )/CTLA.V_2;
    end
end

end