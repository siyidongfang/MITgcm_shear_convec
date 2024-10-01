clear all;close all
fontsize = 20;

fg1=figure(1);
clf;set(gcf,'Color','w','Position',[231 346 600 400])

% tiledlay = tiledlayout(1,2);
nexttile
load('fig4/Ri_topo4.mat')
plot(1./Ri_min,shear_calc_Ri*1e3,'LineWidth',2)
hold on;
load('fig4/Ri_flat.mat')
plot(1./Ri_min,shear_calc_Ri*1e3,'LineWidth',2)
grid on;grid minor;
xlim([0 8.5])
ylim([0 3])

set(gca,'FontSize',fontsize)
legend('Sloping bottom ($\theta=4^\circ$)','Flat bottom ($\theta=0^\circ$)',...
    'Position',[0.4922 0.2860 0.3861 0.13],'interpreter','latex');
legend('boxoff');
ylabel('Tidal shear $\Lambda$ (10$^{-3}\,$s $^{-1}$)','interpreter','latex');
xlabel('${R_i}_\mathrm{min}^{-1}$: Reciprocal of minimum Ri','interpreter','latex');
title('${R_i}_\mathrm{min}^{-1}$ vs Background Shear','interpreter','latex','FontSize',fontsize+3);


% clear Ri_min shear_calc_Ri
% 
% nexttile
% load('fig4/Ri_topo4_harm.mat')
% plot(1./Ri_mean,shear_calc_Ri*1e3,'LineWidth',2)
% hold on;
% load('fig4/Ri_flat_harm.mat')
% plot(1./Ri_mean,shear_calc_Ri*1e3,'LineWidth',2)
% grid on;grid minor;
% xlim([0 4.2])
% ylim([0 3])
% set(gca,'FontSize',fontsize)
% legend('Sloping bottom ($\theta=4^\circ$)','Flat bottom ($\theta=0^\circ$)',...
%     'Position', [0.7 0.2889 0.1727 0.14],'interpreter','latex');
% legend('boxoff');
% ylabel('Tidal shear $\Lambda$ (10$^{-3}\,$s $^{-1}$)','interpreter','latex');
% xlabel('${\overline{R_i}}^{-1}$: Reciprocal of time-averaged Ri','interpreter','latex');
% title('Mean ${R_i}$ vs Background Shear','interpreter','latex','FontSize',fontsize+3);



tiledlay.TileSpacing = 'compact';
tiledlay.Padding = 'compact';
% AddLetters2Plots(fg1,'FontSize',fontsize+5,'FontWeight','normal')

% print('-dpng','-r300',['fig_supp_new/figS_RivsShear.png']);

