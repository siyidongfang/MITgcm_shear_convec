clear all;close all

addpath ../observations/
addpath ../analysis/colormaps/
% load_colors;

fg1 = figure(1);
clf;
set(gcf,'Color','w','Position',[100 153 1200 1000])
tiledlay = tiledlayout(4,1);

%--- MAVS2
load('MAVS2_shear_movmean.mat') 
fontsize = 20;

nexttile;
plot(time/86400,shear_linear)
grid on;grid minor;
set(gca,'Fontsize',fontsize);
ylabel('shear (1/s)','Interpreter','latex')
title('Linear-fit shear at MAVS2','Interpreter','latex')
ylim([-1.7 2.2]/1e3)
xlim([0 90])

nexttile;
plot(time(1:Nend)/86400,abs(shear_fit),'LineWidth',1,'Color','k')
ylabel('shear (1/s)','Interpreter','latex')
grid on;grid minor;
set(gca,'Fontsize',fontsize);
title('Shear magnitude at MAVS2 (sinusoidal curve fit using a 4-tidal-cycle moving window)','Interpreter','latex')


%--- MAVS1
load('MAVS1_shear_movmean.mat') 
fontsize = 20;

nexttile;
plot(time/86400,shear_linear)
grid on;grid minor;
set(gca,'Fontsize',fontsize);
ylabel('shear (1/s)','Interpreter','latex')
title('Linear-fit shear at MAVS1','Interpreter','latex')
ylim([-1.7 2.2]/1e3)
xlim([0 90])

nexttile;
plot(time(1:Nend)/86400,abs(shear_fit),'LineWidth',1,'Color','k')
xlabel('Time (days)','Interpreter','latex')
ylabel('shear (1/s)','Interpreter','latex')
grid on;grid minor;
set(gca,'Fontsize',fontsize);
title('Shear magnitude at MAVS1 (sinusoidal curve fit using a 4-tidal-cycle moving window)','Interpreter','latex')


tiledlay.TileSpacing = 'compact';
tiledlay.Padding = 'compact';

addpath ~/MITgcm_shear_convec/figures/
AddLetters2Plots(fg1,'FontSize',fontsize+5,'FontWeight','normal','Direction','TopDown')

% print('-dpng','-r300','fig_supp_new/figS_obs_shear_movmean.png');
