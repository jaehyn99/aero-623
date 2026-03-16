function [GPhix, GPhiy] = gbasis(xn, xy)
% evaluates basis function gradients at xy

% obtain 1D basis functions and gradients
Phi1x  = basis1d(xn, xy(:,1));
Phi1y  = basis1d(xn, xy(:,2));
GPhi1x = gbasis1d(xn, xy(:,1));
GPhi1y = gbasis1d(xn, xy(:,2));

% return tensor product
nn = length(xn);
GPhix = zeros(size(xy,1), nn*nn);
GPhiy = zeros(size(xy,1), nn*nn);
for q = 1:size(xy,1),
  GPx = GPhi1x(q,:)' *  Phi1y(q,:);  % this is an nn by nn matrix
  GPy =  Phi1x(q,:)' * GPhi1y(q,:);  % this is an nn by nn matrix
  GPhix(q,:) = reshape(GPx, 1, nn*nn);
  GPhiy(q,:) = reshape(GPy, 1, nn*nn);
end

%------------------------------
function GPhi = gbasis1d(xn, x)
% evaluates basis function gradients at x

n = length(x);
order = length(xn)-1;
GPhi = zeros(n, order+1);
for p=0:order,
  B = glagrange(xn, p+1,x);
  GPhi(:,p+1) = B';
end