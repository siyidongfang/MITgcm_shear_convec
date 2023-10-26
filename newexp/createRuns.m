
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

Nx = 200;

Hmax = 1500;
Hsurface = 1000; 
Ntop = 50;
dz_const = 6;
Nr = round((Hmax-Hsurface)/dz_const) + Ntop + 1;

run_type = 'spin'; %%% select from 'init','spin','prod' for initialize run with very small time step, spin-up run for 10 to 20 days, and product run,

%%% Name of the simulation
% exp_name = createRunName (Atide,randtopog_height,randtopog_length,Nr,Nx,run_type)

exp_name = 'f0_noise1e-6'


newexp(batch_name,exp_name,Atide,randtopog_height,randtopog_length,Nr,Nx,run_type)








