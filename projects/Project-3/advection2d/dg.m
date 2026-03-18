function U = dg(L,V,T,N,p,Nt)
% Solves unsteady 2D advection problem using DG in space
%
%    u_t + V dot grad u = 0,    x,y lie in [0,L]^2
%
% L  = length of computational domain; square in x,y
% V  = x,y components of advection velocity
% T  = extent of time domain
% N  = number of elements in each dimension
% p  = order of solution approximation in space
% Nt = number of time steps
%
% Time integration is performed using RK4 with a constant dt.

% build a quadrilateral mesh
[Vert, E2N, I2E, h, normals] = buildmesh(L,N);

% number of elements
nelem = N*N;

% number of coefficients per element (assume tensor-product basis)
nn = (p+1)*(p+1);

% quadrature points sufficient to integrate 2*p in 1D and 2D
nq1 = p+1;                       % number of 1d points
[sq,wq1] = lgwt(nq1,0,1);        % these are 1d points
[nq, xyq, wq] = quad2d(nq1,0,1); % these are 2d points

% nodes for Lagrange basis -- equally-spaced
xn = linspace(0,1,p+1)';
xyn = zeros(nn,2); k = 0;
for j=1:p+1, for i=1:p+1, k = k+1; xyn(k,:) = [xn(i), xn(j)]; end; end

% number of degrees of freedom
Ndof = nelem*nn;

% initialize the solution
U = zeros(Ndof,1);
for j=1:N,
  for i=1:N,
    m = (j-1)*N+i;
    xyglob = h*[(i-1) + xyn(:,1), (j-1) + xyn(:,2)];
    U(nn*(m-1) + [1:nn]) = exactsol(L,xyglob,0);
  end
end


% calculate inverse mass matrix for one element
Phi = basis(xn,xyq);
M   = h*h*Phi'*diag(wq)*Phi;
iM  = inv(M);

% wrap up data for residual evaluation into one structure
resdata.E2N     = E2N;
resdata.I2E     = I2E;
resdata.normals = normals;
resdata.h       = h;
resdata.V       = V;
resdata.xn      = xn;
resdata.xyq     = xyq;
resdata.wq      = wq;
resdata.sq      = sq;
resdata.wq1     = wq1;

% Time stepping
dt = T/Nt;
for it = 1:Nt,
  % RK4
  R = calcresidual(U           , resdata);
  F0 = -massinvmult(iM, R);
  R = calcresidual(U + .5*dt*F0, resdata);
  F1 = -massinvmult(iM, R);
  R = calcresidual(U + .5*dt*F1, resdata);
  F2 = -massinvmult(iM, R);
  R = calcresidual(U +    dt*F2, resdata);
  F3 = -massinvmult(iM, R);
  U = U + (dt/6)*(F0 + 2*F1 + 2*F2 + F3);
end


%----------------------------------------
% evaluates f = M^{-1}*R
function f = massinvmult(iM, R)
nn    = size(iM,1);
nelem = size(R,1)/nn;
f     = zeros(size(R));
for m=1:nelem, % loop over elements
  Im  = (m-1)*nn+[1:nn];
  f(Im) = iM*R(Im);
end

    







