
clear;
addpath ../transition/


fontsize = 18;

fg1 = figure(1);
clf;
set(gcf,'Color','w','Position', [-46 153 1339 800])

tiledlay = tiledlayout(3,3);
x = 0:0.1:5;
xt = 1;
yt = 0.5;
str = '$R=1.7$';
nexttile
plot(x,sin(x))
grid on;grid minor;set(gca,'FontSize',fontsize)
title('Time-dependent frequency $\sigma(\hat t)$','Interpreter','latex','FontSize',fontsize+5);
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+3);

nexttile
plot(x,sin(x+1))
grid on;grid minor;set(gca,'FontSize',fontsize)
title('Horizontal vorticity perturbation $\zeta(\hat t)$','Interpreter','latex','FontSize',fontsize+5);
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+3);

nexttile
plot(x,sin(x+2))
grid on;grid minor;set(gca,'FontSize',fontsize)
title('Turbulent kinetic energy $\sim0.5*\zeta^2(\hat t)$','Interpreter','latex','FontSize',fontsize+5);
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+3);

nexttile
plot(x,sin(x))
grid on;grid minor;set(gca,'FontSize',fontsize)
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+3);

nexttile
plot(x,sin(x+1))
grid on;grid minor;set(gca,'FontSize',fontsize)
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+3);

nexttile
plot(x,sin(x+2))
grid on;grid minor;set(gca,'FontSize',fontsize)
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+3);

nexttile
plot(x,sin(x))
grid on;grid minor;set(gca,'FontSize',fontsize)
xlabel('Time, $\hat t/(2\pi)$ (Tidal cycles)','Interpreter','latex');
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+3);

nexttile
plot(x,sin(x+1))
grid on;grid minor;set(gca,'FontSize',fontsize)
xlabel('Time, $\hat t/(2\pi)$ (Tidal cycles)','Interpreter','latex');
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+3);

nexttile
plot(x,sin(x+2))
grid on;grid minor;set(gca,'FontSize',fontsize)
xlabel('Time, $\hat t/(2\pi)$ (Tidal cycles)','Interpreter','latex');
text(xt,yt,str,'Interpreter','latex','FontSize',fontsize+3);


tiledlay.TileSpacing = 'compact';
tiledlay.Padding = 'compact';




AddLetters2Plots(fg1,'FontSize',fontsize+5,'FontWeight','normal')



% print('-dpng','-r300',['fig_supp_new/figS_ana_flat1.png']);
