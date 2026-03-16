function [nq, xyq, wq] = quad2d(nq1,a,b);
% nq1^2 quadrature points on [a,b]^2

[sq,wq1] = lgwt(nq1,a,b);   % these are 1d points
nq = nq1*nq1;
xyq = zeros(nq,2);
wq  = zeros(nq,1);
k = 0;
for j=1:nq1, for i=1:nq1,
    k = k+1;
    xyq(k,:) = [sq(i), sq(j)];
    wq(k) = wq1(i)*wq1(j);
end; end

