
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

batch_name = 'exps_kv5e-5'


%%% Input parametersd
Atide = 0;

randtopog_height = 0; %%% 10
randtopog_length = 0; %%% 1000


Nx = 1000;
% Nr = 500;
Nr = 800;

% Hmax = 950;
% dz_const = 3;
% Hsurface = 1002; 
% Ntop = 120;
% Nr = round((Hmax-Hsurface)/dz_const) + Ntop;
% Nr = round(Hmax/dz_const);

run_type = 'spin'; %%% select from 'init','spin','prod' for initialize run with very small time step, spin-up run for 10 to 20 days, and product run,

%%% Name of the simulation
% exp_name = createRunName (Atide,randtopog_height,randtopog_length,Nr,Nx,run_type)

%flat [0.1 0.4 0.6 0.8 1.0:0.1:2.2 1.455]*1e-03;
%topo4 [0.1 0.4 0.6 0.8 1.0:0.1:2.0 2.07 1.335]*1e-03;
% topo4 [0.1 0.4:0.2:1.4 1.5:0.1:2.0 2.07 1.335]*1e-03;

Shear = 0.9e-3;
exp_name = ['zeroCenter_topo4_H800Lx3k_s' num2str(Shear*1e3,'%.3f') 'dz1dx3sm100']
newexp(batch_name,exp_name,Atide,randtopog_height,randtopog_length,Nr,Nx,run_type,Shear)



% % % %%% Critical-slope
% % % om = 2*pi/43200;
% % % f0 = 1.18e-4;
% % % N2 = 1e-6;
% % % r_c = sqrt((om^2-f0^2)/(N2-om^2))
% % % slope_c = atand(r_c)
% % % critical_slope = 4.912113091736391






