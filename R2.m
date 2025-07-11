function R = R2()
dev = zeros(4,1);

vehicle_date  = [100; 270.588; 531.373; 713.725; 1133.333]*8e4;
PDL1_data = [100; 217.647; 421.568; 607.843; 964.705]* 8e4;
CTLA4_data  = [100; 227.451; 372.549; 564.706; 1005.882]* 8e4;
combine_data = [100; 203.922; 278.431; 407.843; 539.216]* 8e4;
real = [vehicle_date, PDL1_data, CTLA4_data, combine_data]';

name = {'vehicle','anti-PD-L1', 'anti-CTLA-4', 'combine'};
for i = 1:4
    filename = ['data/basic/',name{i},'.dat'];
    A = load(filename);
    data = [A(1,2), A(351,2), A(701,2), A(1051,2), A(1401,2)];
    dev(i) = calculate_R2(data, real(i,:));
end
R = [dev; sum(dev)]';

end

function result = calculate_R2(simulation, real)

RSS = sum((real - simulation).^2);      % 残差平方和
TSS = sum((real - mean(real)).^2); % 总平方和
result = 1 - (RSS / TSS);

end