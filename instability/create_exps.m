
clear
close all;
addpath ../analysis/colormaps/

FigureIsVisible = 'off';

Shear = 1*0.8e-3; 
Hshear = 300;
N = 1e-3;
topo = 4;
Ptide = 43200;
dz = 0.01;             % dimensionless vertical grid spacing
dt = 0.002;
NTtide = 2;
m1km = 1000;
lambda = 0.5*m1km;


topo_parm = [0 2  6 8 10 12];
N_parm = [0.1 0.5  2  3 4 5 6 7]*1e-3;
Shear_parm = [0 0.2 0.4 0.6   1 1.2 1.4 1.6 1.8 2 2.2 2.4]*1e-3; 
lambda_parm = [0.005 0.01 0.05 0.1 0.5 1 5 10 50 100]*m1km;


for Nexp = 2:length(Shear_parm)
    lambda = lambda_parm(Nexp)
    % N = N_parm(Nexp)
    % topo = topo_parm(Nexp)
    % Shear = Shear_parm(Nexp)
    expdir = ['exps_instability/topo' num2str(topo) '_Pt' num2str(Ptide) ...
        '_N' num2str(N) '_Hs' num2str(Hshear) '_S' num2str(Shear) ...
        '_L' num2str(lambda) '_dz' num2str(dz) '_dt' num2str(dt) ...
        '_NT' num2str(NTtide) '/']
    
    mkdir(expdir);
    
    numerical;
    decompose;
end
