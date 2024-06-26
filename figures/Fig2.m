
clear;close all;
addpath ../analysis/colormaps/
fontsize = 15;
load_colors;

figure(1)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 900 950]);

%%% coordinate
ax1 = subplot('position',[.03 .785 .3 .2]);
annotation('textbox',[0 0.993 0.15 0.01],'String','A','FontSize',fontsize+3,'fontweight','normal','LineStyle','None');
annotation('textbox',[0.07 0.993 0.3 0.01],'String','Slope-aligned coordinate','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
imshow('fig2/coordinate.png')

%--- Load MITgcm simulation
addpath ../analysis/
addpath ../analysis/functions/
expname = 'topo0_H500_s0.0016dz1dx3ln200n-20sm100_kv2e-4';
expdir = '../exps_hires/';
loadexp;
filename = [expdir expname '/RMSE.mat'];
load(filename)
load('fig2.mat')
YLIM = [0 300];

%%% TKE time series
ax2 = subplot('position',[0.435 0.815 0.505 0.16]);
annotation('textbox',[0.38 0.993 0.15 0.01],'String','B','FontSize',fontsize+3,'fontweight','normal','LineStyle','None');
plot(time_h/12,log(pe)/2,'LineWidth',2);
hold on;
plot(time_h/12,log(ke)/2,'LineWidth',2);
plot(xxplot(fit_span)/12, y_fit(fit_span),'k--','LineWidth',2);
xlabel('Time (tidal cycles)')
ylabel('log(energy)')
set(gca,'Fontsize',fontsize)
h2 = legend('Turbulent potential energy','Turbulent kinetic energy','Linear fit',...
    'Fontsize',fontsize+1,'Position',[0.69 0.825 0.2285 0.0661]);
title('Normalized turbulent energy in the shear layer','Fontsize',fontsize+3)
grid on;grid minor;
hold on;
ylim([-38 2])


%%% Velocity
ax3 = subplot('position',[0.07 0.62 0.87 0.125]);
annotation('textbox',[0 0.755 0.15 0.01],'String','C','FontSize',fontsize+3,'fontweight','normal','LineStyle','None');
pcolor(time_tidal,zz-botZ,uu_timeseries');
hold on;shading interp;
contour(time_tidal,zz-botZ,uu_timeseries',[0.15:0.15:0.75],'color',darkgray)
contour(time_tidal,zz-botZ,uu_timeseries',[0 0],'color',darkgray,'LineWidth',1)
contour(time_tidal,zz-botZ,uu_timeseries',[-0.75:0.15:-0.15],'--','color',darkgray)
colormap(redblue);set(gca,'Fontsize',fontsize);
title('$u$','Fontsize',fontsize+4,'interpreter','latex','Position',[15,288])
clim([-0.6 0.6])
ylabel('HAB (m)');
ylim(YLIM)
h3=colorbar(ax3);
set(h3,'Position',[0.95    0.62    0.008    0.11]);
set(get(h3,'Title'),'String',{'$\ \ \ \ (\mathrm{m/s})$'},'interpreter','latex');
set(gca,'xtick',[])
colormap(cmocean('balance'));


%%% Temperature
ax4 = subplot('position',[0.07 0.46 0.87 0.13]);
annotation('textbox',[0 0.595 0.15 0.01],'String','D','FontSize',fontsize+3,'fontweight','normal','LineStyle','None');
pcolor(time_tidal,zz-botZ,tt_timeseries');
hold on;shading interp;
contour(time_tidal,zz-botZ,uu_timeseries',[0.15:0.15:0.75],'color',darkgray)
contour(time_tidal,zz-botZ,uu_timeseries',[0 0],'color',darkgray,'LineWidth',1)
contour(time_tidal,zz-botZ,uu_timeseries',[-0.75:0.15:-0.15],'--','color',darkgray)
title('$\theta^\prime$','Fontsize',fontsize+4,'interpreter','latex','Position',[15,288])
ylabel('HAB (m)');
clim([-0.1 0.1]/5);
set(gca,'Fontsize',fontsize);
ylim(YLIM)
h4=colorbar(ax4);
set(h4,'Position',[0.95    0.46   0.008    0.11]);
set(get(h4,'Title'),'String',{'$\ \ \ \ (^\circ \mathrm{C})$'},'interpreter','latex');
set(gca,'xtick',[])


%%% N2
ax5 = subplot('position',[0.07 0.3 0.87 0.13]);
annotation('textbox',[0 0.435 0.15 0.01],'String','E','FontSize',fontsize+3,'fontweight','normal','LineStyle','None');
pcolor(time_tidal,zz-botZ,(N2_timeseries)')
hold on;
contour(time_tidal,zz-botZ,(N2_timeseries)',[0 0],'Color','c','LineWidth',1);
contour(time_tidal,zz-botZ,uu_timeseries',[0.15:0.15:0.75],'color',darkgray)
contour(time_tidal,zz-botZ,uu_timeseries',[0 0],'color',darkgray,'LineWidth',1)
contour(time_tidal,zz-botZ,uu_timeseries',[-0.75:0.15:-0.15],'--','color',darkgray)
shading interp;
xlabel('Tidal cycles');ylabel('HAB (m)')
set(gca,'Fontsize',fontsize);
title('$N^2$','Fontsize',fontsize+4,'interpreter','latex','Position',[15,285])
clim(([-1 1]+1)/1e6)
ylim(YLIM)
h5=colorbar(ax5);
set(h5,'Position',[0.95    0.3   0.008    0.1]);
set(get(h5,'Title'),'String',{'$\ \ \ \ (1/\mathrm{s}^2)$',''},'interpreter','latex');

%%% Temperature snapshot
ax6 = subplot('position',[0.05 0.05 0.4 0.18]);
annotation('textbox',[0.02 0.2 0.15 0.01],'String','F','FontSize',fontsize+3,'fontweight','normal','LineStyle','None');

%%% N2 snapshot
ax7 = subplot('position',[0.55 0.05 0.4 0.18]);
annotation('textbox',[0.52 0.2 0.15 0.01],'String','G','FontSize',fontsize+3,'fontweight','normal','LineStyle','None');


%%% Save the figure

% print('-djpeg','-r300','fig2/fig2_matlab.png');
