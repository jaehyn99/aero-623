function compare;
% driver routine that compares several DG runs

% problem definitions
L = 1;
V = [1,1];
T = 1;
Nt = 1000; % overkill in most cases

% loop over orders
vN = [4,8,16];
pmax = 3;
E = zeros(pmax+1,length(vN));
for p=2:pmax,
  % loop over numbers of elements
  for k=1:length(vN),
    U = dg(L,V,T, vN(k), p, Nt);
    E(p+1,k) = plotsol(U,L,p);
    title(sprintf('p = %d, N = %d', p, vN(k)), 'fontsize', 24);
    xlabel('x', 'fontsize', 24);
    ylabel('y', 'fontsize', 24);
    fname = sprintf('solN%dp%d.png', vN(k), p);
    print('-dpng',fname);
    %unix(sprintf('epstopdf %s', fname));
  end;
end

% rate plot
ctype = {'ko-', 'bo-', 'ro-', 'go-', 'mo-', 'co-'};
figure(2); clf;
for p=0:pmax,
  dof = vN*(p+1);
  loglog(1./sqrt(dof), E(p+1,:), ctype{p+1}, 'linewidth', 2); hold on;
  rate = log2(E(p+1,end-1)/E(p+1,end));
  slegend{p+1} = sprintf('p = %d, rate = %.2g\n', p, rate);
end
set(gca, 'fontsize', 16);
xlabel('1/sqrt(DOF)', 'fontsize', 24);
ylabel('L_2 solution error norm', 'fontsize', 24);
h = legend(slegend, 4); set(h, 'fontsize', 16);
grid on;
fname = 'rates.png';
print('-dpng',fname);
%unix(sprintf('epstopdf %s', fname));
