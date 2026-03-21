function Phi = basis(xn, xy);
% evaluates basis functions at xy

% obtain 1D basis functions
Phix = basis1d(xn, xy(:,1));
Phiy = basis1d(xn, xy(:,2));

% return tensor product
nn = length(xn);
Phi = zeros(size(xy,1), nn*nn);
for q = 1:size(xy,1),
  P = Phix(q,:)' * Phiy(q,:);  % this is an nn by nn matrix
  Phi(q,:) = reshape(P, 1, nn*nn);
end

