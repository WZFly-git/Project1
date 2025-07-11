function Mat = solve(y0, T)

%% select durg parameter

Cell = y0;
[~, Other] = fun(T(1), Cell);
m = length(T);      n = length([Cell; Other])+1;
Mat = zeros(m,n);
Mat(1,:) = [T(1); Cell; Other]';

%% four order L-K to solve ODE
global md
for k = 2:length(T)
    K1 = fun(T(k-1), Cell);
    K2 = fun(T(k-1) + md.h/2, Cell + K1*md.h/2);
    K3 = fun(T(k-1) + md.h/2, Cell + K2*md.h/2);
    K4 = fun(T(k-1) + md.h, Cell + K3*md.h);
    Cell = Cell + (K1 + 2*K2 + 2*K3 + K4)*md.h/6;
    [~, Other] = fun(T(k), Cell);
    Mat(k,:) = [T(k); Cell; Other]';
end


end


