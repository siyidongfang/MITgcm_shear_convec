
clear;close all;

fontsize = 16;

figure(1)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 1000 350]);

load('fig5/fig5_topo0_Ri_eig.mat')
ax1 = axes('position',[0.075 0.15+0.002 0.39 0.77]);
annotation('textbox',[0.015 0.99+0.002 0.15 0.01],'String','A','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
scatter(omega0_sort/omega,sqrt(epsilon_sort)/omega,20,grow_sort,'o', 'filled');
clim([0 0.04])
% colormap([[1 1 1];cmocean('phase')])
colormap([[1 1 1]*0.97;jet])
set(gca,'fontsize',fontsize)
ylabel('$\sqrt\epsilon_0/\omega = \sqrt {\frac{2\omega_0^3}{\omega}\frac{\Lambda}{\tilde N}}\Big/\omega$','Interpreter','latex','FontSize',fontsize+4)
xlabel('$\omega_0 /\omega = (\tilde N \frac{k_0}{m_0})/\omega$','Interpreter','latex','FontSize',fontsize+4);
title('Growth rate (hour$^{-1}$)','Interpreter','latex','FontSize',fontsize+5)
box on;
grid on;grid minor
xlim([0 2.5])

%%
load('fig5/fig5_topo4_Ri_eig.mat')
ax2 = axes('position',[0.55 0.15+0.002 0.39 0.77]);
annotation('textbox',[0.49 0.99+0.002 0.15 0.01],'String','B','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
scatter(omega0_sort/omega,sqrt(epsilon_sort)/omega,20,grow_sort,'o','filled');
clim([0 0.04])
% ylabel('$\sqrt\epsilon/\omega = \sqrt {2\Lambda\omega_0^3/(N\omega^3)}$','Interpreter','latex')
% xlabel('$\sqrt\delta/\omega = \omega_0 /\omega = N k_x/m_z/\omega$','Interpreter','latex')
set(gca,'fontsize',fontsize)
ylabel('$\sqrt\epsilon_0/\omega = \sqrt {\frac{2\omega_0^3}{\omega}\frac{\Lambda}{\tilde N}}\Big/\omega$','Interpreter','latex','FontSize',fontsize+4)
xlabel('$\omega_0 /\omega = (\tilde N \frac{k_0}{m_0})/\omega$','Interpreter','latex','FontSize',fontsize+4)
title('Growth rate (hour$^{-1}$)','Interpreter','latex','FontSize',fontsize+5)
box on;
grid on;grid minor
h3 = colorbar(ax2,'Position',[0.95 0.1400+0.002 0.0160/2 0.7800]);
xlim([0 2.5])
% % print('-dpng','-r300','fig5/fig5_matlab.png');

print('-r300','-depsc2', 'fig5/fig5_matlab_eig');







%%
figure(2)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 1000 350]);

load('fig5/fig5_topo0_Ri_eig.mat')
ax1 = axes('position',[0.075 0.15+0.002 0.39 0.77]);
annotation('textbox',[0.015 0.99+0.002 0.15 0.01],'String','C','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
% scatter(omega0_sort/omega,sqrt(epsilon_sort)/omega,20,grow_sort,'o', 'filled');
scatter(omega0_sort/omega,Ri_sort,20,grow_sort,'o', 'filled');
clim([0 0.04])
% colormap([[1 1 1];cmocean('phase')])
colormap([[1 1 1]*0.97;jet])
set(gca,'fontsize',fontsize)
set(gca, 'YScale', 'log')
% ylabel('$\sqrt\epsilon_0/\omega = \sqrt {\frac{2\omega_0^3}{\omega}\frac{\Lambda}{\tilde N}}\Big/\omega$','Interpreter','latex','FontSize',fontsize+4)
ylabel('${R_i}_\mathrm{min}$','Interpreter','latex','FontSize',fontsize+4)
xlabel('$\omega_0 /\omega = (\tilde N \frac{k_0}{m_0})/\omega$','Interpreter','latex','FontSize',fontsize+4);
title('Growth rate (hour$^{-1}$)','Interpreter','latex','FontSize',fontsize+5)
box on;
grid on;grid minor
xlim([0 2.5])
ylim([3 2000])


load('fig5/fig5_topo4_Ri_eig.mat')
ax2 = axes('position',[0.55 0.15+0.002 0.39 0.77]);
annotation('textbox',[0.49 0.99+0.002 0.15 0.01],'String','D','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
% scatter(omega0_sort/omega,sqrt(epsilon_sort)/omega,20,grow_sort,'o','filled');
scatter(omega0_sort/omega,Ri_sort,20,grow_sort,'o','filled');
clim([0 0.04])
% ylabel('$\sqrt\epsilon/\omega = \sqrt {2\Lambda\omega_0^3/(N\omega^3)}$','Interpreter','latex')
% xlabel('$\sqrt\delta/\omega = \omega_0 /\omega = N k_x/m_z/\omega$','Interpreter','latex')
set(gca,'fontsize',fontsize)
set(gca, 'YScale', 'log')
% ylabel('$\sqrt\epsilon_0/\omega = \sqrt {\frac{2\omega_0^3}{\omega}\frac{\Lambda}{\tilde N}}\Big/\omega$','Interpreter','latex','FontSize',fontsize+4)
ylabel('${R_i}_\mathrm{min}$','Interpreter','latex','FontSize',fontsize+4)
xlabel('$\omega_0 /\omega = (\tilde N \frac{k_0}{m_0})/\omega$','Interpreter','latex','FontSize',fontsize+4)
title('Growth rate (hour$^{-1}$)','Interpreter','latex','FontSize',fontsize+5)
box on;
grid on;grid minor
h3 = colorbar(ax2,'Position',[0.95 0.1400+0.002 0.0160/2 0.7800]);
xlim([0 2.5])
ylim([3 2000])
% % print('-dpng','-r300','fig5/fig5_matlab2.png');

print('-r300','-depsc2', 'fig5/fig5_matlab2_eig');



