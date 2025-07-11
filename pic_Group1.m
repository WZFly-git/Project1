function pic_Group1()

% Group1();

parameter();
control();

global md 
T = 0:md.h:md.endTime_dose;
n = length(T);
Q1 = zeros(1, n);   Q2 = zeros(1, n);

D1 = md.total_d1;   Plan1 = generate(D1);
D2 = md.total_d2;   Plan2 = generate(D2);

close all

%% PK cover
figure('Position', [100,100,1600,800]);
colors = [0.918, 0.239, 0.2;0.929, 0.502, 0.231;0.459, 0.922, 0.278;0.333, 0.733, 0.922];

for r = 1:length(Plan1)
    md.G1 = Plan1{r};
    md.G2 = Plan2{r};
    subplot(2,2,r);
    for l = 1:n
        Q1(l) = Durg1(T(l));
        Q2(l) = Durg2(T(l));
    end
    plot(T, Q1, 'Color', colors(r,:), 'LineWidth', 1);
    grid on; box on; hold on;
    plot(T, Q2, '--', 'Color', colors(r,:), 'LineWidth', 1);
    ylim([0 1])
    xlabel('Time (day)');
    ylabel('Durg concentration (nmol/L)');
    set(gca,'FontSize',12,'FontName','Arial');
end
% print('Figure/durg_plan/Group1/PK', '-dpng', '-r600');

%% therapeutic effect
figure('Position', [1100,100,600,400]);
target_1 = load('data/Durg_Plan/Group1/target_1.dat');
target_2 = load('data/Durg_Plan/Group1/target_2.dat');
data = [target_1, target_2];
b = bar(data);
grid on; box on; hold on;

xticks(1:4);
xticklabels({'Plan 1', 'Plan 2', 'Plan 3', 'Plan 4'});

set(b(1), 'FaceColor', 'flat', 'CData', colors, 'EdgeColor', 'flat');
set(b(2), 'FaceColor', 'flat', 'CData', colors, 'EdgeColor', 'flat', 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);
ylabel('\psi');

set(gca,'FontSize',12,'FontName','Arial');
print('Figure/durg_plan/Group1/target', '-dpng', '-r600');

end



function Group1()

parameter();
control();

global md 
T = 0:md.h:md.endTime_dose;
y0 = initialization();
D1 = md.total_d1;
D2 = md.total_d2;

%% Q1
choose_par_durg(2);
Mat = solve(y0, T);
C_star = Mat(2101,2);

Plan1 = generate(D1);
n_plan1 = length(Plan1);
target_1 = zeros(n_plan1, 1);
for i = 1:n_plan1
    md.G1 = Plan1{i};
    Mat = solve(y0, T);
    C = Mat(2101,2);
    target_1(i) = (C_star - C) / C_star;
end


parameter();
control();

choose_par_durg(3);
Mat = solve(y0, T);
C_star = Mat(2101,2);

Plan2 = generate(D2);
n_plan2 = length(Plan2);
target_2 = zeros(n_plan2, 1);
for i = 1:n_plan2
    md.G2 = Plan2{i};
    Mat = solve(y0, T);
    C = Mat(2101,2);
    target_2(i) = (C_star - C) / C_star;
end

writematrix(target_1, 'data/Durg_Plan/Group1/target_1.dat');
writematrix(target_2, 'data/Durg_Plan/Group1/target_2.dat');

end



function P = generate(D)

plan1 = {[0,1,2,3,4], [D, D, D, D, D]./5};
plan2 = {[0,2,4,6,8], [D, D, D, D, D]./5};
plan3 = {[0,3,6,9,12], [D, D, D, D, D]./5};
plan4 = {[0,4,8,12,16],[D, D, D, D, D]./5};

P = {plan1,plan2,plan3,plan4};

end


% function P = gen1()
% 
% global md
% 
% D = md.total_d2;
% plan1 = {[0,2,4,6,8,10,12,14], ones(1,8)*D/8};
% plan2 = {[0,3,6,9,12], ones(1,5)*D/5};
% plan3 = {[0,2,4,6,8,10,12,14], [1/2,1/4,1/8,1/16,1/32,1/64,1/128,1/128]*D};
% plan4 = {[0,3,6,9,12], [1/2,1/4,1/8,1/16,1/16]*D};
% plan5 = {[0,1], ones(1,2)*D/2};
% 
% plan8 = {0,D};
% 
% 
% 
% P = {plan1,plan2,plan3,plan4,plan5,plan6,plan7,plan8};
% 
% end


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
