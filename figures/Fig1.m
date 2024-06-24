
% clear;close all;
addpath ../analysis/colormaps/
fontsize = 16;
load_colors;

%%% Data for plotting bathymetry
addpath ../observations/topography/
ncfname = 'blt_canyon_mb_qc.nc';
lat = ncread(ncfname,'lat')';
lon = ncread(ncfname,'lon')';
z = ncread(ncfname,'z');
load Rockall_gebco.mat
% elev(elev>=0)=NaN;

figure(1)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 1400 800]);


%%% Canyon bathymetry
ax1 = subplot('position',[0.06 0.53 0.25 0.25]);
annotation('textbox',[0.025 0.99 0.15 0.01],'String','A','FontSize',fontsize+2,'fontweight','normal','LineStyle','None');
surf(lon,lat,z'/1000,'EdgeColor','None');
shading flat;
set(gca, 'ZDir','reverse')
% colormap(flip(haxby));
% colormap((cmocean('rain')))

clim([min(min(z))/1000-max(max(z))/1000 max(max(z))/1000])
% zlabel('z (km)')
set(gca,'FontSize',fontsize,'TickDir', 'in','TickLength',[0.1, 0.014]);
xlabel('Longitude','Position',[-12.1 54.07 4]);
ylabel('Latitude','Position',[-12.1 54.2286 3.7223])
zlabel('Depth (km)');
title('Canyon','FontSize',fontsize+2,'Position',[-12.2 54.17 -2.1]);
view([-78 58]);
set(get(gca,'xlabel'),'rotation',64);
set(get(gca,'ylabel'),'rotation',-3.3);
% set(get(gca,'xlabel'),'rotation',38);
% set(get(gca,'ylabel'),'rotation',-6.5);
% h1 = colorbar(ax1,'YDir', 'reverse' );
% set(h1,'Position',[0.245 0.625 0.008 0.1]);
% set(get(h1,'Title'),'String','(km)');
axis tight;
ax = gca;
ax.GridColor = [0.1, 0.1, 0.1]*5; 
% zlim([0.45 3])
box off
% set(gca,'view', [-63.1680 63.4425])
xticks([-12.2:0.1:-11.8]);
% yticks([54.1:0.05:54.35])


%%% Rockall Trough
ax11 = subplot('position',[0.06 0.8 0.25 0.17]);
annotation('textbox',[0.025 0.74 0.15 0.01],'String','B','FontSize',fontsize+2,'fontweight','normal','LineStyle','None');

surf(ele_lon,ele_lat,-elev'/1000,'EdgeColor','None');
set(gca, 'ZDir','reverse')
shading flat;
colormap(cmocean('topo','negative'));
clim([-4.5 4.5])
set(gca, 'ZDir','reverse')
% set(gca,'view', [-63.1680 63.4425])
view([-78 58]);
zlabel('(km)')
xlabel('Longitude','Position',[-18.3013 49.1 6]);
ylabel('Latitude','Position',[-22 54.7165 1]);
zlabel('Depth (km)')
set(get(gca,'xlabel'),'rotation',56);
set(get(gca,'ylabel'),'rotation',-3);
axis tight;box off;
set(gca,'FontSize',fontsize,'TickDir', 'in','TickLength',[0.1, 0.014]);
title('Rockall Trough','FontSize',fontsize+2,'Position',[-5 58 -0.4]);
xticks([-19:3:-6]);
annotation('ellipse',[0.225 0.875 0.021 0.01],'Color','red','LineWidth',2,'rotation',70);
annotation('line','LineStyle','--','Color',red1,'LineWidth',1,'Position', [0.06 0.56 0.165 0.322]);% %%% Observed velocity
annotation('line','LineStyle','--','Color',red1,'LineWidth',1,'Position', [0.31 0.75 -0.08 0.161]);% %%% Observed velocity
% ax2 = subplot('position',[0.38 0.555 0.25 0.4]);
% annotation('textbox',[0.375 0.99 0.15 0.01],'String','B','FontSize',fontsize+2,'fontweight','normal','LineStyle','None');

% %%% Linear-fit velocity
% ax3 = subplot('position',[0.73 0.555 0.25 0.4]);
% annotation('textbox',[0.725 0.99 0.15 0.01],'String','C','FontSize',fontsize+2,'fontweight','normal','LineStyle','None');
% 
% %%% Temperature
% ax4 = subplot('position',[0.03 0.055 0.25 0.4]);
% annotation('textbox',[0.025 0.49 0.15 0.01],'String','D','FontSize',fontsize+2,'fontweight','normal','LineStyle','None');
% 
% %%% Reconstructed dbdz using the observed velocity
% ax5 = subplot('position',[0.38 0.055 0.25 0.4]);
% annotation('textbox',[0.375 0.49 0.15 0.01],'String','E','FontSize',fontsize+2,'fontweight','normal','LineStyle','None');
% 
% %%% Reconstructed dbdz using the linear-fit velocity
% ax6 = subplot('position',[0.73 0.055 0.25 0.4]);
% annotation('textbox',[0.725 0.49 0.15 0.01],'String','F','FontSize',fontsize+2,'fontweight','normal','LineStyle','None');
