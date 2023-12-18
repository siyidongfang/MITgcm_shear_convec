
clear
close all;
addpath ../analysis/colormaps/

FigureIsVisible = 'off';

Shear = 1e-3; 
N = 1e-3;
topo = 4;
Ptide = 43200;
dz = 0.002;         % dimensionless vertical grid spacing
dt = 0.002;
m1km = 1000;
Uconst = 0;
Hdepth = 1500;
delta = Hdepth;
NTtide = 5;

dz_parm = [0.3:0.1:2]*0.01;
topo_parm = [1e-20 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
N_parm = [1e-20 0.01 0.05 0.1 0.25 0.5 0.75 1 2 3 4 5 6 7 8 9 10]*1e-3;
Shear_parm = [0 1e-20 0.05 0.1:0.1:1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6]*1e-3; 
lambda_parm = [0.001 0.0025 0.005 0.0075 0.01 0.025 0.05 0.075 0.1 0.25 0.5 0.75 1 2.5 5 7.5 10 25 50 75 100 250 500 1000]*m1km;
% lambda_parm = [0.001 0.0025 0.005 0.0075 0.01 0.025 0.05 0.075 0.1 0.25 0.5 0.75 1 1.2 1.5 1.7 2 2.5 3 3.5 4 5 7.5 10 25 50 75 100 250 500 1000 2000 5000 10000]*m1km;
Ptide_parm = [0.5:0.5:5 10000]*43200;

for ne =7:length(lambda_parm)
% for ne =7:9
% for ne =1:length(Shear_parm)
% for ne =15:21
% for ne = 1:length(dz_parm)
    % dz = dz_parm(ne)
    % dt = dz
    % lambda = 400;
    lambda = lambda_parm(ne)
    kx = 2*pi/(lambda/delta);
    Shear = 0.0013;
    % topo = 1e-100;
    % kx = 0
    % N = N_parm(ne)
    % topo = topo_parm(ne)
    % Shear = Shear_parm(ne)
    % Ptide = Ptide_parm(ne)

    expdir = ['exps_lambda_S0.0013_0.002/H' num2str(Hdepth) '_topo' num2str(topo) '_Pt' num2str(Ptide) ...
        '_N' num2str(N) '_S' num2str(Shear) ...
        '_lambda' num2str(lambda) '_dz' num2str(dz) '_dt' num2str(dt) '_RK4+AB3/']

    mkdir(expdir);

    NOdiffusion = false;  %%% Exclude diffusion/dissipation
    useParallel = false;  %%% Parallel computing
    useRK4 = false;       %%% Use Tourth-order Runge-Kutta method
    useAB3 = true;       %%% Use Third-order Adams-Bashforth method

    numerical;
    decompose;
end
