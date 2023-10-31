
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

Nx = 320;

Hmax = 1500;
dz_const = 4;
% Hsurface = 1000; 
% Ntop = 200;
% Nr = round((Hmax-Hsurface)/dz_const) + Ntop + 1;
Nr = round(Hmax/dz_const)+1;

run_type = 'spin'; %%% select from 'init','spin','prod' for initialize run with very small time step, spin-up run for 10 to 20 days, and product run,

%%% Name of the simulation
% exp_name = createRunName (Atide,randtopog_height,randtopog_length,Nr,Nx,run_type)

exp_name = ['ushear5.5e-4Hs300_Nr',num2str(Nr),'Nx',num2str(Nx) '_' run_type];
% exp_name = 'extgw_nonhydro_vinit_shear200_f0_noise1e-3'
% exp_name = 'test_gw5_nogw'

newexp(batch_name,exp_name,Atide,randtopog_height,randtopog_length,Nr,Nx,run_type)








