%%%
%%% analytical_km.m
%%%
%%% Solve the buoyancy equation analytically after transformation
%%% Assume cross-isobath and slope-normal wavenumbers k and m

clear all;

syms zeta(t) N2 m0 k shear dzeta0 zeta0

Coeff = -N2/(1+(m0/k-shear*t)^2);
ode =  diff(zeta,t,2)+ Coeff*zeta ==0

Dzeta = diff(zeta,t);

cond = [Dzeta(0) == dzeta0; zeta(0) == zeta0];
% ySol(t) = dsolve(ode,cond,'ExpansionPoint',0);
ySol(t) = dsolve(ode,cond);
% ySol = simplify(ySol);

% figure(1)
% fplot(ySol,[0 86400])

% t=tt;
% dd = double(ySol(t));
% plot(tt/43200,dd)

