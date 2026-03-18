function u = exactsol(L, xyglob, t);
% exact solution function

u = exp(-100*((xyglob(:,1)/L-.5).^2 + (xyglob(:,2)/L-.5).^2));