
%%% createRuns.m
%%%
%%% Creates simulations using newexp with specified input parameters.
%%%

addpath /Users/ysi/MITgcm_shear_convec/analysis;
addpath /Users/ysi/MITgcm_shear_convec/analysis/functions/;
addpath /Users/ysi/MITgcm_shear_convec/analysis/colormaps/;
close all;clear;

batch_name = 'exps_test';
% batch_name = 'exps_shear_convec/spin_down';

%%% Input parametersd
Atide = 0;

randtopog_height = 0; %%% 10
randtopog_length = 0; %%% 1000

Nx = 150;

Hmax = 1500;
dz_const = 3;
Hsurface = 1000; 
Ntop = 150;
Nr = round((Hmax-Hsurface)/dz_const) + Ntop + 1;
% Nr = round(Hmax/dz_const)+1;

run_type = 'spin'; %%% select from 'init','spin','prod' for initialize run with very small time step, spin-up run for 10 to 20 days, and product run,

%%% Name of the simulation
% exp_name = createRunName (Atide,randtopog_height,randtopog_length,Nr,Nx,run_type)

exp_name = 's1.6e-3_noise1e-11'
% exp_name = 'test_shear1.5e-3_phase0.5pi_NOnoise'

% exp_name = ['nh_shear9Hs250_Nr',num2str(Nr),'Nx',num2str(Nx) '_' run_type '_restore'];
% exp_name = ['nh_shear9Hs250_Nr',num2str(Nr),'Nx',num2str(Nx) '_' run_type '_noTideInW'];

% exp_name = 'nonhydro_test'

newexp(batch_name,exp_name,Atide,randtopog_height,randtopog_length,Nr,Nx,run_type)



% % % %%% Critical-slope
% % % om = 2*pi/43200;
% % % f0 = 1.18e-4;
% % % N2 = 1e-6;
% % % r_c = sqrt((om^2-f0^2)/(N2-om^2))
% % % slope_c = atand(r_c)
% % % critical_slope = 4.912113091736391






