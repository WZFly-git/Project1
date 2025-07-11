function choose_par_durg(index)

global PDL CTLA 

Para = load('input/key_para.dat');

par_durg1 = PDL.alpha2_default;
par_durg2 = CTLA.K_Q2_default;


if index == 1
    PDL.alpha2 = par_durg1;
    CTLA.K_Q2 = par_durg2;
elseif index == 2
    PDL.alpha2 = Para(9);
    CTLA.K_Q2 = par_durg2;
elseif index == 3
    PDL.alpha2 = par_durg1;
    CTLA.K_Q2 = Para(14);
else
    PDL.alpha2 = Para(9);
    CTLA.K_Q2 = Para(14);
end


end