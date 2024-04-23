%%%
%%% analytical_km.m
%%%
%%% Solve the buoyancy equation analytically after transformation
%%% Assume cross-isobath and slope-normal wavenumbers k and m

%% Small shear limit
% % % % % % % clear all;
% % % % % % % 
% % % % % % % 
% % % % % % % syms zeta(t) dzetadt0 zeta0 a1 a2 omega
% % % % % % % 
% % % % % % % ode =  diff(zeta,t,t)+ a1*(1+2*a2*sin(omega*t)+3*a2^2*sin(omega*t)^2)*zeta ==0
% % % % % % % 
% % % % % % % Dzeta = diff(zeta,t);
% % % % % % % cond = [Dzeta(0) == dzetadt0; zeta(0) == zeta0];
% % % % % % % 
% % % % % % % % ode =  diff(zeta,t,t)+ a1*zeta ==0
% % % % % % % 
% % % % % % % ySol(t) = dsolve(ode);
% % % % % % % % ySol(t) = dsolve(ode,cond);
% % % % % % % % ySol = simplify(ySol);


%% Large shear limit
clear all;


syms zeta(t) dzetadt0 zeta0 a1 omega shear N

a1 = 2*N^2*omega^2/shear^2
ode =  diff(zeta,t,t)+ a1/cos(2*omega*t)*zeta ==0

Dzeta = diff(zeta,t);
cond = [Dzeta(0) == dzetadt0; zeta(0) == zeta0];

ySol(t) = dsolve(ode);
% ySol(t) = dsolve(ode,cond);



