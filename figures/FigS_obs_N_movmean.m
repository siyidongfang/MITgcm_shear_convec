clear all;close all

addpath ../observations/
addpath ../analysis/colormaps/
load_colors;

fontsize = 18;
M2_period =12*3600+25*60; 
YLIM = [-1.5 13]/1e6;

fg1 = figure(1);
clf;
set(gcf,'Color','w','Position',[100 153 850 1000])
tiledlay = tiledlayout(4,1);

%--- MAVS2
load('MAVS2_N2.mat') 
notnan = ~isnan(N2_zavg);
N2_zavg = N2_zavg(notnan);
smooth_N2_zavg = smooth_N2_zavg(notnan);
time = 1:length(N2_zavg);
N2_movmean = movmean(smooth_N2_zavg,4*M2_period);
nexttile;
plot(time/86400,N2_zavg,'Color',gray,'LineWidth',1)
hold on;
plot(time/86400,smooth_N2_zavg,'k--','LineWidth',1.5)
plot(time/86400,N2_movmean,'LineWidth',2)
grid on;grid minor;
set(gca,'Fontsize',fontsize);
ylabel('($10^{-6}$\,s$^{-2}$)','Interpreter','latex')
title('MAVS2: vertical buoyancy gradient $\overline{\partial_{\tilde z}\mathcal{B}}^{\tilde z}$ (bottom 224\,m)','Interpreter','latex')
ylim(YLIM);xlim([0 14.5])
legend('$\overline{\partial_{\tilde z}\mathcal{B}}^{\tilde z}$: using unsmoothed temperature',...
    '$\overline{\partial_{\tilde z}\mathcal{B}}^{\tilde z}$: using smoothed temperature',...
    'Background $\tilde N^2$, moving average with a 4-tidal-cycle sliding window',...
    'Interpreter','latex','Position', [0.23 0.895 0.6412 0.07],'Fontsize',fontsize+1);
legend('boxoff')

nexttile;
clear N2_movmean N2_zavg smooth_N2_zavg time_temp
load('MAVS2_N2_100m.mat') 
notnan = ~isnan(N2_zavg);
N2_zavg = N2_zavg(notnan);
smooth_N2_zavg = smooth_N2_zavg(notnan);
time = 1:length(N2_zavg);
N2_movmean = movmean(smooth_N2_zavg,4*M2_period);
plot(time/86400,N2_zavg,'Color',gray,'LineWidth',1)
hold on;
plot(time/86400,smooth_N2_zavg,'k--','LineWidth',1.5)
plot(time/86400,N2_movmean,'LineWidth',2)
grid on;grid minor;
set(gca,'Fontsize',fontsize);
ylabel('($10^{-6}$\,s$^{-2}$)','Interpreter','latex')
ylim(YLIM);xlim([0 14.5])
title('MAVS2: vertical buoyancy gradient $\overline{\partial_{\tilde z}\mathcal{B}}^{\tilde z}$ (bottom 96\,m)','Interpreter','latex')


%%

%--- MAVS1
nexttile;
clear N2_movmean N2_zavg smooth_N2_zavg time_temp
load('MAVS1_N2.mat') 
time = 1:length(time_temp);
N2_movmean = movmean(smooth_N2_zavg,4*M2_period);
plot(time/86400,N2_zavg,'Color',gray,'LineWidth',1)
hold on;
plot(time/86400,smooth_N2_zavg,'k--','LineWidth',1.5)
plot(time/86400,N2_movmean,'LineWidth',2)
grid on;grid minor;
set(gca,'Fontsize',fontsize);
ylabel('($10^{-6}$\,s$^{-2}$)','Interpreter','latex')
ylim(YLIM);xlim([0 32])
title('MAVS1: vertical buoyancy gradient $\overline{\partial_{\tilde z}\mathcal{B}}^{\tilde z}$ (bottom 224\,m)','Interpreter','latex')


nexttile;
clear N2_movmean N2_zavg smooth_N2_zavg time_temp
load('MAVS1_N2_100m.mat') 
time = 1:length(time_temp);
N2_movmean = movmean(smooth_N2_zavg,4*M2_period);
plot(time/86400,N2_zavg,'Color',gray,'LineWidth',1)
hold on;
plot(time/86400,smooth_N2_zavg,'k--','LineWidth',1.5)
plot(time/86400,N2_movmean,'LineWidth',2)
grid on;grid minor;
set(gca,'Fontsize',fontsize);
ylabel('($10^{-6}$\,s$^{-2}$)','Interpreter','latex')
xlabel('Time (days)','Interpreter','latex')
ylim(YLIM);xlim([0 32])
title('MAVS2: vertical buoyancy gradient $\overline{\partial_{\tilde z}\mathcal{B}}^{\tilde z}$ (bottom 96\,m)','Interpreter','latex')


tiledlay.TileSpacing = 'compact';
tiledlay.Padding = 'compact';

addpath ~/MITgcm_shear_convec/figures/
AddLetters2Plots(fg1,'FontSize',fontsize+5,'FontWeight','normal','Direction','TopDown')

print('-dpng','-r300','fig_supp_new/figS_obs_N_movmean.png');
