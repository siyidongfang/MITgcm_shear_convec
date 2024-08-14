
%%% createRuns.m
%%%
%%% Creates simulations using newexp with specified input parameters.
%%%

addpath /Users/ysi/MITgcm_shear_convec/analysis;
addpath /Users/ysi/MITgcm_shear_convec/analysis/functions/;
addpath /Users/ysi/MITgcm_shear_convec/analysis/colormaps/;
close all;clear;

% batch_name = 'exps_topo4_noCori'
% batch_name = 'exps_flat_noCori'

batch_name = 'exps_topo4_kv2e-5'


%%% Input parametersd
Atide = 0;

randtopog_height = 0; %%% 10
randtopog_length = 0; %%% 1000


Nx = 1000;
Nr = 500;
% Nr = 800;

% Hmax = 950;
% dz_const = 3;
% Hsurface = 1002; 
% Ntop = 120;
% Nr = round((Hmax-Hsurface)/dz_const) + Ntop;
% Nr = round(Hmax/dz_const);

run_type = 'spin'; %%% select from 'init','spin','prod' for initialize run with very small time step, spin-up run for 10 to 20 days, and product run,

%%% Name of the simulation
% exp_name = createRunName (Atide,randtopog_height,randtopog_length,Nr,Nx,run_type)

%[0.1:0.2:0.7 0.8:0.1:2.5]*1e-03;
% [0.3 0.5 0.8 1.0 1.2 1.4]*1e-03;
Shear =  0.8e-3;
% exp_name = ['topo4_H500Lx3k_s' num2str(Shear*1e3) 'dz1dx3n-20sm100_kv2e-5']
exp_name = ['topo4_H500Lx3k_s0.8dz1dx3n-20sm100_kv8e-5']
newexp(batch_name,exp_name,Atide,randtopog_height,randtopog_length,Nr,Nx,run_type,Shear)



% % % %%% Critical-slope
% % % om = 2*pi/43200;
% % % f0 = 1.18e-4;
% % % N2 = 1e-6;
% % % r_c = sqrt((om^2-f0^2)/(N2-om^2))
% % % slope_c = atand(r_c)
% % % critical_slope = 4.912113091736391






