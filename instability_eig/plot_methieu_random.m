
clear;close all;

fontsize = 16;



figure(1)
clf;   
set(gcf,'Color','w');

subplot(2,2,1)
load('../figures/fig5/fig5_topo0_noDiff_eig_new.mat')
grow_sort(grow_sort<=0)=NaN;
scatter(omega0_sort/omega,Ri_sort,20,grow_sort,'o', 'filled');
clim([0 0.4])
colormap(jet)
set(gca,'fontsize',fontsize)
set(gca, 'YScale', 'log')
ylabel('${R_i}_\mathrm{min}$','Interpreter','latex','FontSize',fontsize+4)
xlabel('$\omega_0 /\omega = (\tilde N \frac{k_0}{m_0})/\omega$','Interpreter','latex','FontSize',fontsize+4);
title('Growth rate (hour$^{-1}$)','Interpreter','latex','FontSize',fontsize+4)
box on;
grid on;grid minor
xlim([0 5])
ylim([0.05 1000])
colorbar

subplot(2,2,2)
load('../figures/fig5/fig5_topo0_noDiff_eig_imag.mat')
grow_sort(grow_sort<=0)=NaN;
scatter(omega0_sort/omega,Ri_sort,20,grow_sort,'o', 'filled');
clim([0 0.4])
colormap(jet)
set(gca,'fontsize',fontsize)
set(gca, 'YScale', 'log')
ylabel('${R_i}_\mathrm{min}$','Interpreter','latex','FontSize',fontsize+4)
xlabel('$\omega_0 /\omega = (\tilde N \frac{k_0}{m_0})/\omega$','Interpreter','latex','FontSize',fontsize+4);
title('Growth rate (hour$^{-1}$)','Interpreter','latex','FontSize',fontsize+4)
box on;
grid on;grid minor
xlim([0 5])
ylim([0.05 1000])
colorbar

% subplot(2,2,2)
% load('../figures/fig5/fig5_topo0_noDiff_eig_random.mat')
% grow_sort(grow_sort<=0)=NaN;
% scatter(omega0_sort/omega,Ri_sort,20,grow_sort,'o', 'filled');
% clim([0 0.4])
% colormap(jet)
% set(gca,'fontsize',fontsize)
% set(gca, 'YScale', 'log')
% ylabel('${R_i}_\mathrm{min}$','Interpreter','latex','FontSize',fontsize+4)
% xlabel('$\omega_0 /\omega = (\tilde N \frac{k_0}{m_0})/\omega$','Interpreter','latex','FontSize',fontsize+4);
% title('Growth rate (hour$^{-1}$)','Interpreter','latex','FontSize',fontsize+4)
% box on;
% grid on;grid minor
% xlim([0 5])
% ylim([0.05 1000])
% colorbar
% 
% 
% subplot(2,2,3)
% load('../figures/fig5/fig5_topo0_noDiff_eig_random2.mat')
% grow_sort(grow_sort<=0)=NaN;
% scatter(omega0_sort/omega,Ri_sort,20,grow_sort,'o', 'filled');
% clim([0 0.4])
% colormap(jet)
% set(gca,'fontsize',fontsize)
% set(gca, 'YScale', 'log')
% ylabel('${R_i}_\mathrm{min}$','Interpreter','latex','FontSize',fontsize+4)
% xlabel('$\omega_0 /\omega = (\tilde N \frac{k_0}{m_0})/\omega$','Interpreter','latex','FontSize',fontsize+4);
% title('Growth rate (hour$^{-1}$)','Interpreter','latex','FontSize',fontsize+4)
% box on;
% grid on;grid minor
% xlim([0 5])
% ylim([0.05 1000])
% colorbar
% 
% 
% 
% subplot(2,2,4)
% load('../figures/fig5/fig5_topo0_noDiff_eig_random3.mat')
% grow_sort(grow_sort<=0)=NaN;
% scatter(omega0_sort/omega,Ri_sort,20,grow_sort,'o', 'filled');
% clim([0 0.4])
% colormap(jet)
% set(gca,'fontsize',fontsize)
% set(gca, 'YScale', 'log')
% ylabel('${R_i}_\mathrm{min}$','Interpreter','latex','FontSize',fontsize+4)
% xlabel('$\omega_0 /\omega = (\tilde N \frac{k_0}{m_0})/\omega$','Interpreter','latex','FontSize',fontsize+4);
% title('Growth rate (hour$^{-1}$)','Interpreter','latex','FontSize',fontsize+4)
% box on;
% grid on;grid minor
% xlim([0 5])
% ylim([0.05 1000])
% colorbar
% 
% 
% 
% 
% 
