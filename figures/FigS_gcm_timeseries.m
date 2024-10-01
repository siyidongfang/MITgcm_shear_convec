
clear;close all;
addpath ../analysis/colormaps/
fontsize = 15;
load_colors;

load('fig_supp/FigS_gcm_timeseries.mat')

tidx = tidx/10;
YLIM = [0 400];

figure(1)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 800 900]);

ax1 = subplot('position',[0.07 0.85 0.85 0.12]);
annotation('textbox',[0 0.99 0.15 0.01],'String','A','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(tidx,zz-botZ,tt_timeseries');
hold on;
contour(tidx,zz-botZ,tt_timeseries',[-0.405:0.015:0.405],'Color',darkgray,'LineWidth',0.5);
% contour(tidx,zz-botZ,N2_timeseries',[0 0],'Color','w','LineWidth',0.2);
hold off;
shading flat
colormap(cmocean('balance'))
clim([-0.02 0.28]);
ylim(YLIM)
% freezeColors
set(gca,'Fontsize',fontsize);
xticks([])
ylabel('HAB (m)','Interpreter','latex');
title('Potential temperature $T$','Interpreter','latex','Fontsize',fontsize+4,'Position',[24 395 0]);
h1=colorbar(ax1,'Position',[0.93    0.85+0.01   0.008    0.105]);
set(get(h1,'Title'),'String',{'$\ \ (^\circ\mathrm{C})$'},'interpreter','latex','FontSize',fontsize);
xlim([0 48])
freezeColors(ax1)

ax2 = subplot('position',[0.07 0.7 0.85 0.12]);
annotation('textbox',[0 0.84 0.15 0.01],'String','B','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(tidx,zz-botZ,log10(epsilon_timeseries)');
hold on;shading flat
clim([-11 -8])
ylim(YLIM)
colormap(cmocean('thermal'))
freezeColors(ax2)
set(gca,'Fontsize',fontsize);
xticks([])
ylabel('HAB (m)','Interpreter','latex');
title('The dissipation rate of turbulent kinetic energy $\epsilon=\nu(\nabla \mathbf{u})^2$','Fontsize',fontsize+4,'interpreter','latex','Position',[24 395 0]);
h2=colorbar(ax2,'Position',[0.93    0.7+0.01   0.008    0.105]);
set(get(h2,'Title'),'String',{'$\log(\mathrm{m^2/s^3})$'},'interpreter','latex','FontSize',fontsize);
xlim([0 48])


ax3 = subplot('position',[0.07 0.55 0.85 0.12]);
annotation('textbox',[0 0.69 0.15 0.01],'String','C','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(tidx,zz-botZ,log10(chi_timeseries)');
hold on;shading flat
clim([-11 -9.4])
ylim(YLIM)
colormap(cmocean('haline'))
% colormap hot
freezeColors(ax3);
set(gca,'Fontsize',fontsize);
xticks([])
ylabel('HAB (m)','Interpreter','latex');
title('The dissipation rate of temperature variance $\chi=\kappa(\nabla T)^2$','Fontsize',fontsize+4,'interpreter','latex','Position',[24 395 0]);
h3=colorbar(ax3,'Position',[0.93    0.55+0.01   0.008    0.105]);
set(get(h3,'Title'),'String',{'$\log(\mathrm{^\circ C^2/s})$'},'interpreter','latex','FontSize',fontsize);
xlim([0 48])


Nplot = N2_timeseries;
Nplot(N2_timeseries<0)=0;
Nplot=log10(sqrt(abs(Nplot)));
ax4 = subplot('position',[0.07 0.4 0.85 0.12]);
annotation('textbox',[0 0.54 0.15 0.01],'String','D','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(tidx,zz-botZ,Nplot');
shading flat
hold on;
% contour(tidx,zz-botZ,N2_timeseries',[0 0],'Color','c','LineWidth',0.2);
hold off;
clim([-3.2 -2.8])
% clim([-2 2]*1e-6);
ylim(YLIM)
colormap(cmocean('haline'))
freezeColors(ax4)
set(gca,'Fontsize',fontsize);
xticks([])
ylabel('HAB (m)','Interpreter','latex');
title('Absolute buoyancy frequency $\mathcal{N} = \vert\partial_{\tilde z}\mathcal{B}\vert ^{1/2}$','Fontsize',fontsize+4,'interpreter','latex','Position',[24 395 0]);
h4=colorbar(ax4,'Position',[0.93    0.4+0.01   0.008    0.105]);
set(get(h4,'Title'),'String',{'$\log(\mathrm{1/s})$'},'interpreter','latex','FontSize',fontsize);
xlim([0 48])


ax5 = subplot('position',[0.07 0.25 0.85 0.12]);
annotation('textbox',[0 0.39 0.15 0.01],'String','E','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(tidx,zz-botZ,uu_timeseries'*cosd(4)-ww_timeseries'*sind(4));
hold on;shading flat
ylim(YLIM)
clim([-0.399 0.399])
colormap(cmocean('balance'))
set(gca,'Fontsize',fontsize);
xticks([])
freezeColors(ax5)
ylabel('HAB (m)','Interpreter','latex');
title('Horizontal velocity $\tilde u= u\cos\theta-w\sin\theta$','Fontsize',fontsize+4,'interpreter','latex','Position',[24 390 0]);
h5=colorbar(ax5,'Position',[0.93    0.25+0.01   0.008    0.105]);
set(get(h5,'Title'),'String',{'$\ \ \ \ (\mathrm{m/s})$'},'interpreter','latex','FontSize',fontsize);
xlim([0 48])


ax6 = subplot('position',[0.07 0.1 0.85 0.12]);
annotation('textbox',[0 0.24 0.15 0.01],'String','F','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(tidx,zz-botZ,ww_timeseries'*cosd(4)+uu_timeseries'*sind(4));
hold on;shading flat
ylim(YLIM)
clim([-0.075 0.075])
colormap(cmocean('balance'))
set(gca,'Fontsize',fontsize);
freezeColors(ax6)
ylabel('HAB (m)','Interpreter','latex');
xlabel('Time (hours)','Interpreter','latex','Fontsize',fontsize+3);
t1 = title('Vertical velocity $\tilde w= u\sin\theta+w\cos\theta$','Fontsize',fontsize+4,'interpreter','latex','Position',[24 390 0]);
h6=colorbar(ax6,'Position',[0.93    0.1+0.01   0.008    0.105]);
set(get(h6,'Title'),'String',{'$\ \ \ \ (\mathrm{m/s})$'},'interpreter','latex','FontSize',fontsize);
xlim([0 48])


print('-dpng','-r300',['fig_supp/figS_gcm_timeseries_matlab.png']);



% figure(2)
% clf;   
% set(gcf,'Color','w');
% scrsz = get(0,'ScreenSize');
% set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 800 900]);
% 
% ax2 = subplot('position',[0.07 0.7 0.85 0.12]);
% h2=colorbar(ax2,'Position',[0.93    0.7+0.01   0.008    0.105]);
% set(gca,'Fontsize',fontsize);
% xticks([])
% clim([-11 -8])
% colormap(cmocean('thermal'))
% 
% print('-dpng','-r300',['fig_supp/figS_gcm_timeseries_colorbar2.png']);



% figure(3)
% clf;   
% set(gcf,'Color','w');
% scrsz = get(0,'ScreenSize');
% set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 800 900]);
% 
% 
% ax3 = subplot('position',[0.07 0.55 0.85 0.12]);
% h3=colorbar(ax3,'Position',[0.93    0.55+0.01   0.008    0.105]);
% colormap(cmocean('haline'))
% set(gca,'Fontsize',fontsize);
% xticks([])
% clim([-11 -9.4])
% 
% print('-dpng','-r300',['fig_supp/figS_gcm_timeseries_colorbar3.png']);
% 
% 
% 
% figure(4)
% clf;   
% set(gcf,'Color','w');
% scrsz = get(0,'ScreenSize');
% set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 800 900]);
% ax4 = subplot('position',[0.07 0.4 0.85 0.12]);
% colormap(cmocean('haline'))
% set(gca,'Fontsize',fontsize);
% h4=colorbar(ax4,'Position',[0.93    0.4+0.01   0.008    0.105]);
% clim([-3.2 -2.8])
% 
% print('-dpng','-r300',['fig_supp/figS_gcm_timeseries_colorbar4.png']);
