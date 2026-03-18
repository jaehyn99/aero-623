function compare;
% driver routine that compares several DG runs

% problem definitions
L = 1;
a = 1;
T = 1;
Nt = 1000; % overkill in most cases

ctype = {'b-', 'r-', 'c-', 'm-'};   % colors for plotting
xg = linspace(0,L,500);            % global x for plotting exact

% loop over orders
vN = [8,16,32,64];
pmax = 4;
E = zeros(pmax+1,length(vN));
for p=0:pmax,
  figure(p+1); clf;
  % loop over numbers of elements
  for k=1:length(vN),
    U = dg(L,a,T, vN(k), p, Nt);
    E(p+1,k) = plotsol(U,L,a,p, ctype{k}); hold on;
    slegend{k} = sprintf('N = %d', vN(k));
  end;
  plot(xg, exactsol(L,a,xg,0), 'k-', 'linewidth', 2);
  title(sprintf('p = %d', p), 'fontsize', 24);
  xlabel('x', 'fontsize', 24);
  ylabel('u', 'fontsize', 24);
  grid on; set(gca, 'fontsize', 16);
  h = legend(slegend, 'location','northeast'); set(h, 'fontsize', 16);
  fname = sprintf('solp%d.eps', p);
  print('-depsc2',fname);
  unix(sprintf('epstopdf %s', fname));
end

% rate plot
ctype = {'ko-', 'bo-', 'ro-', 'go-', 'mo-', 'co-'};
figure(pmax+2); clf;
for p=0:pmax,
  dof = vN*(p+1);
  loglog(1./dof, E(p+1,:), ctype{p+1}, 'linewidth', 2); hold on;
  rate = log2(E(p+1,end-1)/E(p+1,end));
  slegend{p+1} = sprintf('p = %d, rate = %.2g\n', p, rate);
end
set(gca, 'fontsize', 16);
xlabel('1/DOF', 'fontsize', 24);
ylabel('L_2 solution error norm', 'fontsize', 24);
h = legend(slegend, 'location','southeast'); set(h, 'fontsize', 16);
grid on;
fname = 'rates.eps';
print('-depsc2',fname);
unix(sprintf('epstopdf %s', fname));
