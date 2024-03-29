%%%
%%% analytical_km.m
%%%
%%% Solve the buoyancy equation analytically after transformation
%%% Assume cross-isobath and slope-normal wavenumbers k and m

clear all;


%%%%% All variables are dimensional variables
N = 1e-3;
N2 = N^2;

shear = 1.3e-3;
Ptide = 43200;
omega = 2*pi/Ptide;
rw = 0.1; %%% ratio of the wavenumbers kx/mz
rs = shear/omega; %%% shear over omega
b0 = 1e-20;
db0 = 1e-25;

NTtide = 20;
dt = 60;
Lt = NTtide*43200; 
Nt = Lt/dt;
tt = dt:dt:Nt*dt;
topo = 0.07;


% syms b(t) db0 b0 rw rs omega N2 topo
syms b(t)

Q = -(rw^2+(rw*rs*sin(omega*t)-1)^2)/((rw*cos(topo)-sin(topo))*N2);
P = rw*cos(topo)-sin(topo) + rw*rs*sin(omega*t)*sin(topo);
dQdt = diff(Q,t);

ode =  Q*diff(b,t,2)+ dQdt*diff(b,t)-P*b ==0

Db = diff(b,t);
cond = [Db(0) == db0; b(0) == b0];
% cond = [Db(0) == 0; b(0) == b0];
% cond = [Db(0) == 0; b(0) == 0];
ySol(t) = dsolve(ode,cond,'ExpansionPoint',0);
% ySol = simplify(ySol);

% figure(1)
% fplot(ySol,[0 86400])

t=tt;
dd = double(ySol(t));
plot(tt/43200,dd)

