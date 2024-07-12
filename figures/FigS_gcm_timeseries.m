
clear;close all;
addpath ../analysis/colormaps/
fontsize = 15;
% load_colors;

load('fig_supp/FigS_gcm_timeseries.mat')
figure(1)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 800 900]);

ax1 = subplot('position',[0.07 0.85 0.87 0.12]);
annotation('textbox',[0 0.99 0.15 0.01],'String','a','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(tidx,zz-botZ,tt_timeseries');
hold on;shading flat
colormap(cmocean('balance'))
colorbar;
clim([0 0.35]);

ax2 = subplot('position',[0.07 0.7 0.87 0.12]);
annotation('textbox',[0 0.84 0.15 0.01],'String','b','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(tidx,zz-botZ,N2_timeseries');
hold on;
contour(tidx,zz-botZ,N2_timeseries',[0 0],'Color','c','LineWidth',0.5);
hold off;
shading flat
colorbar;
clim([0 2]*1e-6);

ax3 = subplot('position',[0.07 0.55 0.87 0.12]);
annotation('textbox',[0 0.69 0.15 0.01],'String','c','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(tidx,zz-botZ,log10(chi_timeseries)');
hold on;shading flat
colorbar;
clim([-12.5 -8.6])

ax4 = subplot('position',[0.07 0.4 0.87 0.12]);
annotation('textbox',[0 0.54 0.15 0.01],'String','d','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(tidx,zz-botZ,log10(epsilon_timeseries)');
hold on;shading flat
colorbar;
clim([-11 -8])

ax5 = subplot('position',[0.07 0.25 0.87 0.12]);
annotation('textbox',[0 0.39 0.15 0.01],'String','e','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(tidx,zz-botZ,uu_timeseries');
hold on;shading flat
colorbar;

ax6 = subplot('position',[0.07 0.05 0.87 0.12]);
annotation('textbox',[0 0.24 0.15 0.01],'String','f','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(tidx,zz-botZ,ww_timeseries');
hold on;shading flat
colorbar;
