function  [K, Other] = fun(t, Cell)

global par Tumor PDL CTLA 

C = Cell(1); D = Cell(2); Tc = Cell(3); Th = Cell(4); Tr = Cell(5);

%% cytokine and path_par updata
I_2 = (par.delta_I2Th)/(par.d_I2)* Th;
I_10 = (par.delta_I10Tr)/(par.d_I10)* Tr + (par.delta_I10C)/(par.d_I10)* C;
I_12 = (par.delta_I12D)/(par.d_I12)* D;
T_beta = (par.delta_TbetaTr)/(par.d_Tbeta)* Tr + (par.delta_TbetaC)/(par.d_Tbeta)* C ;
I_gamma = (par.delta_IgammaTh)/(par.d_Igamma)* Th + (par.delta_IgammaTc)/(par.d_Igamma)* Tc;

P = PDL.rho_P* (Th + Tc);
L = PDL.rho_L* (Th + Tc + par.epsilon_C* C);
Q_1 = Durg1(t);

W = CTLA.rho_W* (Th + Tc);
C_4 = CTLA.rho_C4* (Th + Tc + par.epsilon_r* Tr);
B = CTLA.rho_B* D;
Q_2 = Durg2(t);

%% cell updata
% par_choose();  % 更新关键全局变量的信息
K = zeros(5,1);

K(1) = Tumor.beta_C* (1 - C/Tumor.cell_CM)* C - Tumor.eta_Th* Th* C - Tumor.eta_Tc* Tc* C;

K(2) = par.lambda_DC* ( C/(par.K_C + C) )* par.D0 - par.d_D* D;

M_D = (D^par.Hill_n)/(par.K_D^par.Hill_n + D^par.Hill_n);
M_I12 = (I_12)/(par.K_I12 + I_12);
M_I2 = (I_2)/(par.K_I2 + I_2);

M_TcI10 = (par.K_TcI10)/(par.K_TcI10 + I_10);
M_TcTbeta  = (par.K_TcTbeta)/(par.K_TcTbeta + T_beta);
K(3) = (par.kappa_Tc* M_D* par.TN8* M_I12* M_TcI10 *M_TcTbeta + par.beta_TcI2* M_I2* Tc  )* inhibition(P,L,Q_1)* promotion(W,B,C_4,Q_2) ...
        - par.d_Tc* Tc;

M_ThI10 = (par.K_ThI10)/(par.K_ThI10 + I_10);
M_ThTbeta  = (par.K_ThTbeta)/(par.K_ThTbeta + T_beta);
K(4) = (par.kappa_Th* M_D* par.TN4* M_I12* M_ThI10* M_ThTbeta + par.beta_ThI2* M_I2* Th  )* inhibition(P,L,Q_1)* promotion(W,B,C_4,Q_2)...
        - par.d_Th* Th;

M_Tbeta = (T_beta)/(par.K_Tbeta + T_beta);
M_Igamma = (par.K_TrIgamma)/(par.K_TrIgamma + I_gamma);
K(5) = par.kappa_Tr* M_D* par.TN4* M_Tbeta* M_Igamma - par.d_Tr* Tr;

Other = [I_2; I_10; I_12; T_beta; I_gamma; P; L; Q_1; W; C_4; B; Q_2];

end

%%  inhibition and promotion function
function rate = inhibition(P,L,Q_1)
global PDL
P_L = (PDL.alpha1* P* L)/(1 + PDL.alpha1* P + PDL.alpha2* Q_1);
rate = PDL.K_PL/(PDL.K_PL + P_L);
end

function rate = promotion(W,B,C_4,Q_2)
global CTLA
M_Q = CTLA.K_Q2/(CTLA.K_Q2 + Q_2);
M_C4 = CTLA.K_C4/(CTLA.K_C4 + C_4* M_Q);
rate = (W* B* M_C4)/(CTLA.K_WB + W* B* M_C4);
end

%% Durg function 
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


