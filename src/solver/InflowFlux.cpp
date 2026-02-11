function [Fb] = inflow_flux(UI, Tt, pt, n, alpha)
% calculates boundary flux given interior state (UI),
% total temperature (Tt), total pressure (pt), outward-pointing
% edge normal (n), and flow angle (alpha) [rad].
% gas properties
gam = 1.4; % ratio of specific heats
R = 1; % gas constant

% ensure a normal vector for n
n = n/norm(n);

% inflow direction
nin = [cos(alpha), sin(alpha)];
dn = dot(nin,n);

% properties of interior state
rI = UI(1); % interior density
VI = UI(2:3)/rI; % interior velocity
unI = dot(VI,n); % interior normal speed
pI = (gam-1)*(UI(4) - 0.5*rI*dot(VI,VI)); % interior pressure
cI = sqrt(gam*pI/rI); % interior speed of sound

% Riemann invariant from interior
Jp = unI + 2*cI/(gam-1);

% solve quadratic for boundary Mach number, Mb
A = gam*R*Tt*dn^2 - (gam-1)/2*Jp^2;
B = 4*gam*R*Tt*dn/(gam-1);
C = 4*gam*R*Tt/(gam-1)^2 - Jp^2;
disc = sqrt(B^2-4*A*C);
M = [-B-disc, -B+disc]/(2*A);
Mb = min(M(find(M>0))); % use smallest positive root

% calculate boundary state
Tb = Tt/(1+0.5*(gam-1)*Mb^2); % static temperature
pb = pt*(Tb/Tt)^(gam/(gam-1)); % static pressure
rb = pb/(R*Tb); % static density
cb = sqrt(gam*pb/rb); % speed of sound
Vb = Mb*cb*nin; % velocity
unB = dot(Vb,n); % normal speed
rEb = pb/(gam-1) + 0.5*rb*dot(Vb,Vb); % total energy
% calculate boundary flux
Fb = zeros(4,1);
Fb(1) = rb*unB;
Fb(2:3) = rb*unB*Vb + pb*n;
Fb(4) = (rEb+pb)*unB;