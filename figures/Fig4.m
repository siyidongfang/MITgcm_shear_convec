
clear;close all;

fontsize = 16;

figure(1)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 1000 350]);


ax1 = axes('position',[0.08 0.14 0.39 0.78]);
annotation('textbox',[0.02 0.9 0.15 0.01],'String','b','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');

load('fig4/fig4_topo0.mat')

scatter(omega0_sort/omega,epsilon_sort,20,grow_sort,'o', 'filled');
clim([0 0.3]/10)
colormap([[1 1 1];cmocean('phase')])
ylabel('$\sqrt\epsilon/\omega = \sqrt {2\Lambda\omega_0^3/(N\omega^3)}$','Interpreter','latex')
xlabel('$\sqrt\delta/\omega = \omega_0 /\omega = N k_x/m_z/\omega$','Interpreter','latex')
set(gca,'fontsize',fontsize)
title('Growth rate (1/hour)','Interpreter','latex','FontSize',fontsize+4)
% colorbar;
box on;
grid on;grid minor



load('fig4/fig4_topo4.mat')