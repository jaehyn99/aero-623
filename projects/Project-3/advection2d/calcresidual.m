function R = calcresidual(U, resdata)
% evaluates residual of weak form

% pull off variables from resdata structure
E2N     = resdata.E2N;
I2E     = resdata.I2E;
normals = resdata.normals;
h       = resdata.h;
V       = resdata.V;
xn      = resdata.xn;
xyq     = resdata.xyq;
wq      = resdata.wq;
sq      = resdata.sq;
wq1     = resdata.wq1;

nelem = size(E2N,1);
nn = length(sq)^2;

% create and zero out residual vector
R = zeros(size(U));

% Contributions from element interiors
Phi = basis(xn,xyq);               % basis fcns at quad points
[GPhiX, GPhiY] = gbasis(xn,xyq);   % grad of basis fcns at quad points
GPhix = GPhiX/h; GPhiy = GPhiY/h;  % grad of basis fcns in glob space
for m=1:nelem,
  Im = (m-1)*nn+[1:nn];      % unknown indices on elem m
  Um = U(Im);                % unknown coefficients on elem m
  u = Phi*Um;                % state at quad points
  Fx = u*V(1); Fy = u*V(2);  % flux vector at quad points  
  R(Im) = R(Im) - h*h*(GPhix'*diag(wq)*Fx + GPhiy'*diag(wq)*Fy);
end

% pre-compute basis fcns on each local edge/orientation
EdgePhi = cell(4,2); 
for e=1:4,
  xy = RefEdge2Elem(e, sq);      % counter-clockwise orientation
  EdgePhi{e,1} = basis(xn,xy);
  xy = RefEdge2Elem(e, 1.0-sq);  % clockwise orientation
  EdgePhi{e,2} = basis(xn,xy);
end

% Contributions from interior edges
for ie = 1:size(I2E,1),
  elemL = I2E(ie,1);  elemR = I2E(ie,2); % left/right elements
  edgeL = I2E(ie,3);  edgeR = I2E(ie,4); % left/right loc edges
  n  = normals(ie,:);                    % normal from L to R
  IL = (elemL-1)*nn+[1:nn];  % unkown indices on elemL
  IR = (elemR-1)*nn+[1:nn];  % unkown indices on elemR
  PhiL = EdgePhi{edgeL, 1};  % L basis functions at edge quad points
  PhiR = EdgePhi{edgeR, 2};  % R basis functions at edge quad points
  UL = U(IL);                % unknown coefficients on elemL
  UR = U(IR);                % unknown coefficients on elemR
  uL = PhiL*UL;              % L state at edge quad points
  uR = PhiR*UR;              % R state at edge quad points
  Fhat = fluxfunction(V, uL,uR,n); % upwind flux at edge quad points
  R(IL) = R(IL) + h*PhiL'*diag(wq1)*Fhat;
  R(IR) = R(IR) - h*PhiR'*diag(wq1)*Fhat;
end


%--------------------------------------
function Fhat = fluxfunction(V, uL, uR, n)
% computes upwind flux given left/right states and normal (L->R)
Vdotn = dot(V,n);
if (Vdotn > 0)
  Fhat = uL*Vdotn;
else
  Fhat = uR*Vdotn;
end


%----------------------------------
function xy = RefEdge2Elem(edge, s)
% computes 2D ref-elem coords for 1d edge coords in s
xy = zeros(length(s),2);
switch edge
  case 1
    xy(:,1) = s; xy(:,2) = 0;
  case 2
    xy(:,1) = 1; xy(:,2) = s;
  case 3
    xy(:,1) = 1-s; xy(:,2) = 1;
  case 4
    xy(:,1) = 0; xy(:,2) = 1-s;
end
