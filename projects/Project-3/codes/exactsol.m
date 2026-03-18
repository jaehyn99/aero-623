function u = exactsol(L, a, x, t);
% exact solution function

u = exp(-100*(x/L-.5).^2);