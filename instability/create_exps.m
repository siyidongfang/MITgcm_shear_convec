
clear
close all;
addpath ../analysis/colormaps/

FigureIsVisible = 'on';

Shear = 1e-3; 
N = 1e-3;
topo = 4;
Ptide = 43200;
dz = 0.006;             % dimensionless vertical grid spacing
dt = 0.001;
m1km = 1000;

Hdepth = 1500;
delta = Hdepth;

topo_parm = [1e-20 2 4 6 8 10 12 14 16];
N_parm = [0.01 0.05 0.1 0.25 0.5 0.75 1 2 3 4 5 6 7 8 9 10]*1e-3;
Shear_parm = [1e-20 0.05 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2 2.2 2.4]*1e-3; 
lambda_parm = [0.001 0.0025 0.005 0.0075 0.01 0.025 0.05 0.075 0.1 0.25 0.5 0.75 1 2.5 5 7.5 10 25 50 75 100 250 500 1000 1e7 1e47 1e97]*m1km;

% for Nexp =22:length(lambda_parm)
for Nexp = 4:6
    % lambda = 750;
    lambda = lambda_parm(Nexp)
    kx = 2*pi/(lambda/delta);
    % kx = 0
    % N = N_parm(Nexp)
    % topo = topo_parm(Nexp)
    % Shear = Shear_parm(Nexp)


    expdir = ['exps/H' num2str(Hdepth) '_topo' num2str(topo) '_Pt' num2str(Ptide) ...
        '_N' num2str(N) '_S' num2str(Shear) ...
        '_lambda' num2str(lambda) '_dz' num2str(dz) '_dt' num2str(dt) '/']

    mkdir(expdir);

    %%% Parallel computing
    numerical;
    decompose;
end
