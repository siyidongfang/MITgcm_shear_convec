%%%
%%% analytical_km.m
%%%
%%% Solve the buoyancy equation analytically after transformation
%%% Assume cross-isobath and slope-normal wavenumbers k and m

clear all;

%%%%% All variables are dimensional variables
% N = 1e-3;
% shear = 1.3e-3;
% Ptide = 43200;
% omega = 2*pi/Ptide;
% rw = 1; %%% ratio of the wavenumbers kx/mz
% rs = shear/omega; %%% shear over omega
% b0 = 1e-20;
% db0 = 1e-40;

% NTtide = 5;
% dt = 60;
% Lt = NTtide*43200; 
% Nt = Lt/dt;
% tt = dt:dt:Nt*dt;

% beta1 = 1 + rw^2 + 2*rw*rs*sin(omega*t) + rw^2*rs^2*(sin(omega*t))^2;
% beta2 = 2*rw*shear*cos(omega*t)+2*rw^2*shear*rs*sin(omega*t)*cos(omega*t);
% beta3 = rw^2*N^2;

syms b(t) db0 b0 rw rs omega shear N

% ode = diff((1 + rw^2 + 2*rw*rs*sin(omega*t) + rw^2*rs^2*(sin(omega*t))^2)...
%     *diff(b,t),t) + rw^2*N^2*b== 0;

ode = (1 + rw^2 + 2*rw*rs*sin(omega*t) + rw^2*rs^2*(sin(omega*t))^2) * diff(b,t,2) ...
      + (2*rw*shear*cos(omega*t)+2*rw^2*shear*rs*sin(omega*t)*cos(omega*t)) * diff(b,t) ...
      + rw^2*N^2 * b == 0;

Db = diff(b,t);
cond = [Db(0) == db0; b(0) == b0];
ySol(t) = dsolve(ode,cond,'IgnoreAnalyticConstraints',false)
% ySol = simplify(ySol);

% figure(1)
% fplot(ySol,[0 86400])
