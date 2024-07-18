
clear;close all;

% load_colors;
fontsize = 16;
addpath ../instability/products/
load('GrowthRate_exps_linear_dz0.5.mat')


grow_400m = GrowthRate_Floquet(16,:);
figure(1)
plot(shear_Floquet,grow_400m)

% figure(1)
% plot(shear_Floquet,growth_Floquet)


% load('GrowthRate_exps_tanh_ZeroCenter_dz2.mat')
% figure(1)
% plot(shear_Floquet,growth_Floquet)


% plot
%%% 8 panels 

%--- Flat bottom: growth rate as a function of inverse Ri and log10(lambda_x)

%--- Sloping bottom: growth rate as a function of inverse Ri and log10(lambda_x)

%--- Flat bottom: maximum growth rate as a function of inverse Ri (MITgcm vs theory)

%--- Sloping bottom: maximum growth rate as a function of inverse Ri (MITgcm vs theory)

%--- Flat bottom: buoyancy budget

%--- Sloping bottom: buoyancy budget

%--- Flat bottom: vorticity budget

%--- Sloping bottom: vorticity budget


% print('-dpng','-r200',['fig_supp/FigS_disc_growth.png']);




