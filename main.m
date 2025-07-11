function  main()

%% define global para, initial value and time
parameter();
control();

global md
T = 0:md.h:md.endTime;
y0 = initialization();


%% calculation and save data
for index = 1:4
    choose_par_durg(index); 
    Mat = solve(y0, T);
    save_data(Mat, index);
end


end

function save_data(M, index)

name = {'vehicle','anti-PD-L1', 'anti-CTLA-4', 'combine'};

filename = ['data/basic/',name{index},'.dat'];
fp = fopen(filename, 'w');
[m,n] = size(M);
for i = 1:m
    for j = 1:n
        fprintf(fp, '%8.8e ',M(i,j));
    end
    fprintf(fp,'\n ');
end

fclose all;
end