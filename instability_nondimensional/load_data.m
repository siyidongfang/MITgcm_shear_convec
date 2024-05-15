
N = 1e-3;
topo = 4;
Ptide = 43200;
dz = 0.002;         % dimensionless vertical grid spacing
dt = 0.002;
m1km = 1000;
Uconst = 0; %%% A time-invariant background flow adding to the tidal oscillation
Hdepth = 1500;
delta = Hdepth;
NTtide = 10;

Shear_parm = [0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4]*1e-3; 

expfolder = 'exps_inviscid';

%%% Load data
lambda = 400;
kx = 2*pi/(lambda/delta);

Shear = Shear_parm(ne)
expdir = [expfolder '/H' num2str(Hdepth) '_topo' num2str(topo) '_Pt' num2str(Ptide) ...
    '_N' num2str(N) '_S' num2str(Shear) ...
    '_lambda' num2str(lambda) '_dz' num2str(dz) '_dt' num2str(dt) '/']
load([expdir 'output.mat'])

