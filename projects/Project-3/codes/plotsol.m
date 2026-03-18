function EL2 = plotsol(U,L,a,p, ctype);
% plots a 1D DG solution
% basis is assumed Lagrange on equally-spaced nodes
% output is L2 error norm

nn = p+1;          % number of coefficients per element
N = size(U,1)/nn;  % number of elements

% quadrature points and weights
[xiq,wq] = lgwt(2*nn+3,0,1);

% Lagrange basis nodes
xn = linspace(0,1,nn)';

% plot solution at points xp
xp = linspace(0,1,4*(p+1))';
Phi = basis(xn,xp);  % basis functions at plotting nodes
Phiq = basis(xn,xiq); % basis functions at quadrature points
dx = L/N;
EL2 = 0;
vxg = []; vu = [];
for m=1:N,
  xg = (m-1)*dx + xp*dx;
  % interpolate u
  Um = U((m-1)*nn+[1:nn]);
  u = Phi*Um;
  % store solution for plotting
  vxg = [vxg;xg]; vu = [vu;u];
  % error integration
  xq = (m-1)*dx + xiq*dx;
  uq = Phiq*Um;
  EL2 = EL2 + dx*wq'*(uq-exactsol(L,a,xq,0)).^2;
end
EL2 = sqrt(EL2/L);
plot(vxg, vu, ctype, 'linewidth', 2);


