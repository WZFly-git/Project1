function control()

global md 
md.h = 0.01;
md.endTime = 15;   %21天
md.endTime_dose = 21;

md.total_d1 = 6.7;   md.total_d2 = 6.75; 
md.G1{1} = [0, 3.5, 7, 10.5, 14];      md.G1{2} = [1, 1, 1, 1, 1]* md.total_d1/5;
md.G2{1} = [0, 3.5, 7, 10.5, 14];      md.G2{2} = [1, 1, 1, 1, 1]* md.total_d2/5;



end




