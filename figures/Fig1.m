
% clear;close all;
addpath ../analysis/colormaps/
fontsize = 16;


%%% Data for plotting bathymetry
addpath ../observations/topography/
ncfname = 'blt_canyon_mb_qc.nc';
lat = ncread(ncfname,'lat')';
lon = ncread(ncfname,'lon')';
z = ncread(ncfname,'z');


figure(1)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 1400 800]);

%%% Canyon bathymetry
ax1 = subplot('position',[0.06 0.555 0.25 0.4]);
annotation('textbox',[0.025 0.99 0.15 0.01],'String','A','FontSize',fontsize+2,'fontweight','normal','LineStyle','None');
surf(lon,lat,z'/1000,'EdgeColor','None');
shading flat;set(gca, 'ZDir','reverse')
% colormap(flip(haxby));
colormap((cmocean('rain')))
xlabel('Longitude');
ylabel('Latitude');
zlabel('Depth (km)')
set(gca,'FontSize',fontsize)
title('Canyon topography','FontSize',fontsize+3);
view([-78 58]);
set(get(gca,'xlabel'),'rotation',76);
set(get(gca,'ylabel'),'rotation',-4);
h1 = colorbar(ax1,'YDir', 'reverse' );
set(h1,'Position',[0.245 0.625 0.008 0.1]);
set(get(h1,'Title'),'String','(km)');
axis tight;
ax = gca;
ax.GridColor = [0.1, 0.1, 0.1]*5; 
% grid(gca,'minor') 
zlim([0.45 3])
box on


%%% Observed velocity
ax2 = subplot('position',[0.38 0.555 0.25 0.4]);
annotation('textbox',[0.375 0.99 0.15 0.01],'String','B','FontSize',fontsize+2,'fontweight','normal','LineStyle','None');

%%% Linear-fit velocity
ax3 = subplot('position',[0.73 0.555 0.25 0.4]);
annotation('textbox',[0.725 0.99 0.15 0.01],'String','C','FontSize',fontsize+2,'fontweight','normal','LineStyle','None');

%%% Temperature
ax4 = subplot('position',[0.03 0.055 0.25 0.4]);
annotation('textbox',[0.025 0.49 0.15 0.01],'String','D','FontSize',fontsize+2,'fontweight','normal','LineStyle','None');

%%% Reconstructed dbdz using the observed velocity
ax5 = subplot('position',[0.38 0.055 0.25 0.4]);
annotation('textbox',[0.375 0.49 0.15 0.01],'String','E','FontSize',fontsize+2,'fontweight','normal','LineStyle','None');

%%% Reconstructed dbdz using the linear-fit velocity
ax6 = subplot('position',[0.73 0.055 0.25 0.4]);
annotation('textbox',[0.725 0.49 0.15 0.01],'String','F','FontSize',fontsize+2,'fontweight','normal','LineStyle','None');
