
clear
close all;
addpath ../analysis/colormaps/

FigureIsVisible = 'off';

Shear = 1e-3; 
N = 1e-3;
topo = 4;
Ptide = 43200;
dz = 0.01;         % dimensionless vertical grid spacing
dt = 0.01;
m1km = 1000;
Uconst = 0; %%% A time-invariant background flow adding to the tidal oscillation
Hdepth = 250;
delta = Hdepth;
NTtide = 10;

topo_parm = [1e-20 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
N_parm = [1e-20 0.01 0.05 0.1 0.25 0.5 0.75 1 2 3 4 5 6 7 8 9 10]*1e-3;
% Shear_parm = ([2.3 0.1:0.1:2.2])*1e-3;
% Shear_parm = ([0.1:0.2:2.1])*1e-3;
% Shear_parm = ([0.1:0.4:2.1])*1e-3;
Shear_parm = ([0.2:0.4:2.0])*1e-3;
% lambda_parm = [5 10 50:50:350 400:50:1000 1200:200:5000 6000:1000:20000 30000:10000:100000];
lambda_parm = round(10.^[1.7:0.1:3.4 3.6 3.8 4]/10)*10;
%lambda_parm = round(10.^[1.7:-0.1:1]);
Ptide_parm = [0.5:0.5:5 10000]*43200;

exppath = 'exps_hydro_Nr100/';

% for Nexp_lambda =1:length(lambda_parm)
for Nexp_lambda =4

    lambda = lambda_parm(Nexp_lambda);
    expfolder = [exppath 'lambda' num2str(lambda) '/']

    % for Nexp_shear =1:length(Shear_parm)
    for Nexp_shear =5:6

        Shear = Shear_parm(Nexp_shear)
        expdir = [expfolder 'H' num2str(Hdepth) '_topo' num2str(topo) '_Pt' num2str(Ptide) ...
            '_N' num2str(N) '_S' num2str(Shear) ...
            '_lambda' num2str(lambda) '/'];
        mkdir(expdir);
        
        NOdiffusion = false;  %%% Exclude diffusion/dissipation
        useRK4 = true;      %%% Use Tourth-order Runge-Kutta method
        useAB3 = false;       %%% Use Third-order Adams-Bashforth method

        noBQ2 = false;
        noBQ3 = false;
        noBQ4 = false;
        noZQ2 = false;
        noZQ3 = false;

        hydrostatic = false;

        run_numerical;

    end
    
end
