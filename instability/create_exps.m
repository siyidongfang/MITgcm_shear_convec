
clear
close all;
addpath ../analysis/colormaps/

FigureIsVisible = 'off';

Shear = 1e-3; 
N = 1e-3;
topo = 4;
Ptide = 43200;
dz = 0.005;         % dimensionless vertical grid spacing
dt = 0.005;
m1km = 1000;
Uconst = 0; %%% A time-invariant background flow adding to the tidal oscillation
Hdepth = 300;
delta = Hdepth;
NTtide = 20;

topo_parm = [1e-20 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
N_parm = [1e-20 0.01 0.05 0.1 0.25 0.5 0.75 1 2 3 4 5 6 7 8 9 10]*1e-3;
Shear_parm = [0.1:0.1:2.5]*1e-3; 
lambda_parm = [400:50:1000 1200:200:5000 6000:1000:20000 30000:10000:100000];
Ptide_parm = [0.5:0.5:5 10000]*43200;

expfolder = 'experiments/';

for Nexp_lambda = 1

    lambda = lambda_parm(Nexp_lambda);
    expfolder = [expfolder 'lambda' num2str(lambda) '/']

    for Nexp_shear =1:length(Shear_parm)

        Shear = Shear_parm(Nexp_shear)
        expdir = [expfolder 'H' num2str(Hdepth) '_topo' num2str(topo) '_Pt' num2str(Ptide) ...
            '_N' num2str(N) '_S' num2str(Shear) ...
            '_lambda' num2str(lambda) '/'];
        mkdir(expdir);
        
        NOdiffusion = false;  %%% Exclude diffusion/dissipation
        useRK4 = false;      %%% Use Tourth-order Runge-Kutta method
        useAB3 = true;       %%% Use Third-order Adams-Bashforth method

        numerical;

    end
    
end
