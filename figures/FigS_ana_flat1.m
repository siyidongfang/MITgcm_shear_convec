
clear;
addpath ../transition/


fontsize = 18;

fg1 = figure(1);
clf;
set(gcf,'Color','w','Position', [-46 153 1339 800])

tiledlay = tiledlayout(3,3);
x = 0:0.1:5;

nexttile
plot(x,sin(x),'LineWidth',2)
grid on;grid minor;set(gca,'FontSize',fontsize)
title('Horizontal vorticity perturbation $\zeta(\hat t)$','Interpreter','latex','FontSize',fontsize+5);

nexttile
plot(x,sin(x+1))
grid on;grid minor;set(gca,'FontSize',fontsize)
title('Horizontal vorticity perturbation $\zeta(\hat t)$','Interpreter','latex','FontSize',fontsize+5);

nexttile
plot(x,sin(x+2))
grid on;grid minor;set(gca,'FontSize',fontsize)
title('Turbulent kinetic energy $\sim0.5*\zeta^2(\hat t)$','Interpreter','latex','FontSize',fontsize+5);

nexttile
plot(x,sin(x))
grid on;grid minor;set(gca,'FontSize',fontsize)

nexttile
plot(x,sin(x+1))
grid on;grid minor;set(gca,'FontSize',fontsize)

nexttile
plot(x,sin(x+2))
grid on;grid minor;set(gca,'FontSize',fontsize)

nexttile
plot(x,sin(x))
grid on;grid minor;set(gca,'FontSize',fontsize)
xlabel('Time, $\hat t/(2\pi)$ (Tidal cycles)','Interpreter','latex');

nexttile
plot(x,sin(x+1))
grid on;grid minor;set(gca,'FontSize',fontsize)
xlabel('Time, $\hat t/(2\pi)$ (Tidal cycles)','Interpreter','latex');

nexttile
plot(x,sin(x+2))
grid on;grid minor;set(gca,'FontSize',fontsize)
xlabel('Time, $\hat t/(2\pi)$ (Tidal cycles)','Interpreter','latex');


tiledlay.TileSpacing = 'compact';
tiledlay.Padding = 'compact';




AddLetters2Plots(fg1,'FontSize',fontsize+5,'FontWeight','normal')


% print('-dpng','-r300',['fig_supp_new/figS_ana_flat1.png']);
