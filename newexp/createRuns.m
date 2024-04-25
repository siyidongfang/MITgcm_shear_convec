
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
% Hsurface = 1002; 
% Ntop = 120;
% Nr = round((Hmax-Hsurface)/dz_const) + Ntop;
Nr = round(Hmax/dz_const);

run_type = 'spin'; %%% select from 'init','spin','prod' for initialize run with very small time step, spin-up run for 10 to 20 days, and product run,

%%% Name of the simulation
% exp_name = createRunName (Atide,randtopog_height,randtopog_length,Nr,Nx,run_type)

Shear = 1e-3
exp_name = ['s' num2str(Shear) '_dz3dx20Lx10km_n1e-20Ln400_smoothShear'];
% exp_name = ['s' num2str(Shear) '_test3']

newexp(batch_name,exp_name,Atide,randtopog_height,randtopog_length,Nr,Nx,run_type,Shear)



% % % %%% Critical-slope
% % % om = 2*pi/43200;
% % % f0 = 1.18e-4;
% % % N2 = 1e-6;
% % % r_c = sqrt((om^2-f0^2)/(N2-om^2))
% % % slope_c = atand(r_c)
% % % critical_slope = 4.912113091736391






