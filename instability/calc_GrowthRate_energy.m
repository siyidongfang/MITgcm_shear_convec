

close all;clear;

expdir = '/Volumes/MIT/MITgcm_shear_convec/instability/experiments/lambda';
lambda_parm = [1000 1200 1400 1600 10000 11000 12000 13000 14000 15000 16000 100000];
Shear_parm = [0.1:0.1:0.8]*1e-3;

lambda = lambda_parm(1);
shear = Shear_parm(1); 
expname = ['H300_topo4_Pt43200_N0.001_S' num2str(shear) '_lambda' num2str(lambda) '/'];

load([expdir num2str(lambda) '/' expname '/output.mat'],...
    're_buoy','re_zeta','re_psi','omega','dt','Ptide','Nt','Nr','zz','Nshear',...
    'uuu','www','re_buoyd')

TKE = 0.5*(uuu.^2+0.5*(www(:,1:Nr)+www(:,2:Nr+1)).^2);
TPE = re_buoyd;
