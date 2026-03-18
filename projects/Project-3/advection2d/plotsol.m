function EL2 = plotsol(U,L,p);
% plots a 2D DG solution
% basis is assumed Lagrange on Gauss-Legendre nodes
% output is L2 error norm

nn = (p+1)^2;       % number of coefficients per element
NN = size(U,1)/nn;  % number of elements
N = sqrt(NN);

% quadrature points
[nq, xyq, wq] = quad2d(p+1,0,1);

% nodes for Lagrange basis -- equally-spaced
xn = linspace(0,1,p+1)';
xyn = zeros(nn,2); k = 0;
for j=1:p+1, for i=1:p+1, k = k+1; xyn(k,:) = [xn(i), xn(j)]; end; end

% plot solution at points xyp
s = linspace(0,1,4*(p+1));
np = length(s);
k = 0;
for j=1:np, for i=1:np,
  k = k+1;
  xyp(k,:) = [s(i), s(j)];
end; end;

% patches for plotting
npatch = (np-1)^2;
pv = zeros(npatch,4);
k = 0;
for j=1:np-1, for i=1:np-1,
    k = k+1;
    ii = (j-1)*np + i;
    pv(k,:) = [ii, ii+1, ii+np+1, ii+np];
end; end;


figure(1); clf;
Phi = basis(xn,xyp);
h = L/N;
EL2 = 0;
vxg = []; vu = [];
m = 0;
for mj = 1:N, for mi = 1:N,
  xyg = [(mi-1)*h + xyp(:,1)*h, (mj-1)*h + xyp(:,2)*h];
  % interpolate u
  m = m+1;
  Um = U((m-1)*nn+[1:nn]);
  u = Phi*Um;
  % plot solution on element m
  for l = 1:npatch, 
    patch(xyg(pv(l,:),1), xyg(pv(l,:),2), u(pv(l,:)), 'linestyle', 'none');
  end
  % error integration
  xyng = [(mi-1)*h + xyn(:,1)*h, (mj-1)*h + xyn(:,2)*h];
  EL2 = EL2 + h*h*wq'*(Um-exactsol(L,xyng,0)).^2;
end; end;
EL2 = sqrt(EL2);
axis square;


