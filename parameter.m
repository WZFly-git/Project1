function parameter()

global PDL CTLA par 
Para = load('input/key_para.dat');
select_para(Para);
%% critical durg para 
PDL.F_1 = 0.8;           PDL.V_1 = 2;             PDL.k_e1 = 1;
CTLA.F_2 = 0.8;          CTLA.V_2 = 2;            CTLA.k_e2 = 0.9;

PDL.alpha2_default = 0;         CTLA.K_Q2_default = 1e5;
%% Half-saturation constant
par.K_C = 1e7;          par.K_D = 8e7;

par.K_I2 = 3e5;
par.K_I12 = 500;
par.K_Tbeta = 0.05;

par.K_TcI10 = 0.9;
par.K_TcTbeta = 0.03;
par.K_ThI10 = 1.1;
par.K_ThTbeta = 0.05;
par.K_TrIgamma = 200;

%% relative fix para

% Cell 
par.D0 = 1.94e7;
par.TN8 = 9.39e9;
par.TN4 = 1.88e9;

% growth rate
par.beta_TcI2 = 3.5;    par.beta_ThI2 = 5;

% activation rate
par.lambda_DC = 3;

% apoptosis rate
par.d_D = 0.1;      par.d_Tc = 0.18;
par.d_Th = 0.197;   par.d_Tr = 0.1;

% differentiation rate
par.kappa_Tc = 27.96;
par.kappa_Th = 24.9;
par.kappa_Tr = 1.5;

% Cytokine production rate and degradtion rate
par.delta_I2Th = 1.15e-4;                                       par.d_I2 = 2.376;
par.delta_I10Tr = 1.4e-8;       par.delta_I10C =1.3e-10;        par.d_I10 = 20;
par.delta_I12D = 9e-7;                                          par.d_I12 = 0.672;
par.delta_TbetaTr = 9.24e-10;   par.delta_TbetaC = 1.10e-7;     par.d_Tbeta = 198;
par.delta_IgammaTh = 6.5e-8;    par.delta_IgammaTc = 2.5e-7;    par.d_Igamma = 10;

% Hill para
par.Hill_n = 4;    

% adjustment para
par.epsilon_C = 50; 
par.epsilon_r = 50;

end

function select_para(Para)
% change critical par

global Tumor PDL CTLA
Tumor.beta_C = Para(1);          Tumor.cell_CM = Para(2);     
Tumor.eta_Tc = Para(3)*1e-11;    Tumor.eta_Th = Para(4)*1e-12;      
PDL.rho_P = Para(5);             PDL.rho_L = Para(6);            PDL.K_PL = Para(7);          
PDL.alpha1 = Para(8);            PDL.alpha2 = Para(9);
CTLA.rho_B = Para(10);           CTLA.rho_C4 = Para(11);         CTLA.rho_W = Para(12);
CTLA.K_WB = Para(13);            CTLA.K_Q2 = Para(14);           CTLA.K_C4 = Para(15);  

end