clear; close all;

addpath ../instability_eig/products/

omega = 2*pi/43200;

fontsize = 16;
fg1=figure(1);
clf;set(gcf,'Color','w','Position', [100 153 1000 1000])
tiledlay = tiledlayout(4,2);

nexttile;
load('methieu_M2_topo0_N1.mat','omega0_sort','Ri_sort','grow_sort')
scatter(omega0_sort/omega,Ri_sort,15,grow_sort,'o', 'filled');
clim([0 0.04])
colormap([[1 1 1]*0.97;jet])
set(gca,'fontsize',fontsize)
set(gca, 'YScale', 'log')
ylabel('${R_i}_\mathrm{min}$','Interpreter','latex','FontSize',fontsize+4)
% xlabel('$\omega_0 /\omega = (\tilde N \frac{k_0}{m_0})/\omega$','Interpreter','latex','FontSize',fontsize+4);
title('Growth rate (hour$^{-1}$)','Interpreter','latex','FontSize',fontsize+5)
box on;
grid on;grid minor
xlim([0 2.5])
ylim([3 1000])



nexttile;
load('methieu_M2_topo2_N1.mat','Ri_sort','grow_sort')
scatter(omega0_sort/omega,Ri_sort,15,grow_sort,'o', 'filled');
clim([0 0.04])
colormap([[1 1 1]*0.97;jet])
set(gca,'fontsize',fontsize)
set(gca, 'YScale', 'log')
ylabel('${R_i}_\mathrm{min}$','Interpreter','latex','FontSize',fontsize+4)
% xlabel('$\omega_0 /\omega = (\tilde N \frac{k_0}{m_0})/\omega$','Interpreter','latex','FontSize',fontsize+4);
title('Growth rate (hour$^{-1}$)','Interpreter','latex','FontSize',fontsize+5)
box on;
grid on;grid minor
xlim([0 2.5])
ylim([3 1000])

nexttile;
load('methieu_M2_topo4_N1.mat','Ri_sort','grow_sort')
scatter(omega0_sort/omega,Ri_sort,15,grow_sort,'o', 'filled');
clim([0 0.04])
colormap([[1 1 1]*0.97;jet])
set(gca,'fontsize',fontsize)
set(gca, 'YScale', 'log')
ylabel('${R_i}_\mathrm{min}$','Interpreter','latex','FontSize',fontsize+4)
% xlabel('$\omega_0 /\omega = (\tilde N \frac{k_0}{m_0})/\omega$','Interpreter','latex','FontSize',fontsize+4);
% title('Growth rate (hour$^{-1}$)','Interpreter','latex','FontSize',fontsize+5)
box on;
grid on;grid minor
xlim([0 2.5])
ylim([3 1000])


nexttile;
load('methieu_M2_topo6_N1.mat','Ri_sort','grow_sort')
scatter(omega0_sort/omega,Ri_sort,15,grow_sort,'o', 'filled');
clim([0 0.04])
colormap([[1 1 1]*0.97;jet])
set(gca,'fontsize',fontsize)
set(gca, 'YScale', 'log')
ylabel('${R_i}_\mathrm{min}$','Interpreter','latex','FontSize',fontsize+4)
% xlabel('$\omega_0 /\omega = (\tilde N \frac{k_0}{m_0})/\omega$','Interpreter','latex','FontSize',fontsize+4);
% title('Growth rate (hour$^{-1}$)','Interpreter','latex','FontSize',fontsize+5)
box on;
grid on;grid minor
xlim([0 2.5])
ylim([3 1000])


nexttile;
load('methieu_M2_topo8_N1.mat','Ri_sort','grow_sort')
scatter(omega0_sort/omega,Ri_sort,15,grow_sort,'o', 'filled');
clim([0 0.04])
colormap([[1 1 1]*0.97;jet])
set(gca,'fontsize',fontsize)
set(gca, 'YScale', 'log')
ylabel('${R_i}_\mathrm{min}$','Interpreter','latex','FontSize',fontsize+4)
% xlabel('$\omega_0 /\omega = (\tilde N \frac{k_0}{m_0})/\omega$','Interpreter','latex','FontSize',fontsize+4);
% title('Growth rate (hour$^{-1}$)','Interpreter','latex','FontSize',fontsize+5)
box on;
grid on;grid minor
xlim([0 2.5])
ylim([3 1000])


nexttile;
load('methieu_M2_topo10_N1.mat','Ri_sort','grow_sort')
scatter(omega0_sort/omega,Ri_sort,15,grow_sort,'o', 'filled');
clim([0 0.04])
colormap([[1 1 1]*0.97;jet])
set(gca,'fontsize',fontsize)
set(gca, 'YScale', 'log')
ylabel('${R_i}_\mathrm{min}$','Interpreter','latex','FontSize',fontsize+4)
% xlabel('$\omega_0 /\omega = (\tilde N \frac{k_0}{m_0})/\omega$','Interpreter','latex','FontSize',fontsize+4);
% title('Growth rate (hour$^{-1}$)','Interpreter','latex','FontSize',fontsize+5)
box on;
grid on;grid minor
xlim([0 2.5])
ylim([3 1000])

nexttile;
load('methieu_M2_topo12_N1.mat','Ri_sort','grow_sort')
scatter(omega0_sort/omega,Ri_sort,15,grow_sort,'o', 'filled');
clim([0 0.04])
colormap([[1 1 1]*0.97;jet])
set(gca,'fontsize',fontsize)
set(gca, 'YScale', 'log')
ylabel('${R_i}_\mathrm{min}$','Interpreter','latex','FontSize',fontsize+4)
xlabel('$\omega_0 /\omega = (\tilde N \frac{k_0}{m_0})/\omega$','Interpreter','latex','FontSize',fontsize+4);
% % title('Growth rate (hour$^{-1}$)','Interpreter','latex','FontSize',fontsize+5)
box on;
grid on;grid minor
xlim([0 2.5])
ylim([3 1000])

nexttile;
load('methieu_M2_topo14_N1.mat','Ri_sort','grow_sort')
scatter(omega0_sort/omega,Ri_sort,15,grow_sort,'o', 'filled');
clim([0 0.04])
colormap([[1 1 1]*0.97;jet])
set(gca,'fontsize',fontsize)
set(gca, 'YScale', 'log')
ylabel('${R_i}_\mathrm{min}$','Interpreter','latex','FontSize',fontsize+4)
xlabel('$\omega_0 /\omega = (\tilde N \frac{k_0}{m_0})/\omega$','Interpreter','latex','FontSize',fontsize+4);
% title('Growth rate (hour$^{-1}$)','Interpreter','latex','FontSize',fontsize+5)
box on;
grid on;grid minor
xlim([0 2.5])
ylim([3 1000])


tiledlay.TileSpacing = 'compact';
tiledlay.Padding = 'compact';
AddLetters2Plots(fg1,'FontSize',fontsize+5,'FontWeight','normal')


print('-dpng','-r300','fig_supp_new/figS_sens_methieu_matlab.png');
