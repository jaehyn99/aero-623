function U = dg(L,a,T,N,p,Nt);
% Solves unsteady 1D advection problem using DG in space
%
%    u_t + a*u_x = 0,    x lies in [0,L]
%
% L  = length of computational domain
% a  = advection velocity (> 0)
% T  = extent of time domain
% N  = number of elements
% p  = order of solution approximation in space
% Nt = number of time steps
%
% Time integration is performed using RK4 with a constant dt.

% element size
dx = L/N;

% number of coefficients per element
nn = p+1;

% quadrature points sufficient to integrate 2*p on [0,1]
[xq,wq] = lgwt(p+1,0,1);

% nodes for Lagrange basis
xn = linspace(0,1,nn)';

% number of degrees of freedom
Ndof = N*nn;

% initialize the solution
U = zeros(Ndof,1);
for m=1:N,
  xglob = dx*(m-1) + xn*dx;
  U(nn*(m-1) + [1:nn]) = exactsol(L,a,xglob,0);
end

% calculate inverse mass matrix for one element
Phi = basis(xn,xq);
M   = dx*Phi'*diag(wq)*Phi;
iM  = inv(M);

% calculate A0
GPhi = gbasis(xn,xq);
A0   = -a*GPhi'*diag(wq)*Phi + a*basis(xn,1)'*basis(xn,1);

% calculate AL
AL   = -a*basis(xn,0)'*basis(xn,1);

% Time stepping
dt = T/Nt;
for it = 1:Nt,
  % RK4
  t = (it-1)*dt;  % time at start of stage
  F0 = F(U           , t      , iM,A0,AL);
  F1 = F(U + .5*dt*F0, t+.5*dt, iM,A0,AL);
  F2 = F(U + .5*dt*F1, t+.5*dt, iM,A0,AL);
  F3 = F(U +    dt*F2, t+   dt, iM,A0,AL);
  U = U + (dt/6)*(F0 + 2*F1 + 2*F2 + F3);
end


%----------------------------------------
% evaluates f = -M^{-1}*A*U
function f = F(U, t, iM, A0, AL);
nn = size(iM,1);
N = size(U,1)/nn;
f = zeros(N*nn,1);
for m=1:N, % loop over elements
  mL = m-1; % element to left
  if (mL==0), mL=N; end;  % account for periodic BCs
  ImL = (mL-1)*nn+[1:nn];
  Im  = (m -1)*nn+[1:nn];
  f(Im) = -iM*(AL*U(ImL) + A0*U(Im));
end







