
clear;

load('om_alpha.mat')

Ntide = 30;
dt = 2*pi/300;
tt_hat = 0:dt:2*pi*Ntide;
Nt = length(tt_hat);

%%% Initial condition

%%% Integrate the system forward in time using RK4

%%% Compute the growth rate
