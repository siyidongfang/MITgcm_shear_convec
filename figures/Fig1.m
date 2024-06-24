
% clear;close all;
addpath ../analysis/colormaps/
addpath freezeColors/
fontsize = 16;
load_colors;


figure(1)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 1400 800]);

%%%------ Data for plotting bathymetry
addpath ../observations/topography/
ncfname = 'blt_canyon_mb_qc.nc';
lat = ncread(ncfname,'lat')';
lon = ncread(ncfname,'lon')';
z = ncread(ncfname,'z');
load Rockall_gebco.mat

%%% Canyon bathymetry
ax2 = subplot('position',[0.06 0.53 0.25 0.25]);
annotation('textbox',[0.025 0.99 0.15 0.01],'String','A','FontSize',fontsize+2,'fontweight','normal','LineStyle','None');
surf(lon,lat,z'/1000,'EdgeColor','None');
shading flat;
set(gca, 'ZDir','reverse')
clim([min(min(z))/1000-max(max(z))/1000 max(max(z))/1000])
set(gca,'FontSize',fontsize,'TickDir', 'in','TickLength',[0.1, 0.014]);
xlabel('Longitude','Position',[-12.1 54.065 4]);
ylabel('Latitude','Position',[-12.1 54.2286 3.7223])
zlabel('Depth (km)');
title('Canyon','FontSize',fontsize+2,'Position',[-12.2 54.17 -2.1]);
view([-78 58]);
set(get(gca,'xlabel'),'rotation',64);
set(get(gca,'ylabel'),'rotation',-3.3);
axis tight;
ax = gca;
ax.GridColor = [0.1, 0.1, 0.1]*5; 
box off
xticks([-12.2:0.1:-11.8]);
% ylim([54.1 54.38])
% colormap([[0.1 0.1 0.1]*2;cmocean('rain')]);
colormap(cmocean('tarn'));
% clim([0 4])
% colormap(cmocean('topo','negative'));
clim([-3.5 3.5])
freezeColors;


%%% Rockall Trough
ax1 = subplot('position',[0.06 0.8 0.25 0.17]);
annotation('textbox',[0.025 0.74 0.15 0.01],'String','B','FontSize',fontsize+2,'fontweight','normal','LineStyle','None');
surf(ele_lon,ele_lat,-elev'/1000,'EdgeColor','None');
hold on;
contour3(ele_lon,ele_lat,-elev'/1000,[0 0],'EdgeColor',black);
% contour3(ele_lon,ele_lat,-elev'/1000,[0.45 0.45],'-.','EdgeColor',black,'ShowText','off');
hold off;
set(gca, 'ZDir','reverse')
shading flat;
% colormap(flip(cmocean('topo')));
% colormap(cmocean('topo','negative'));
clim([-3.5 3.5])
% colormap(flip(haxby))
% colormap([[0.1 0.1 0.1]*2;cmocean('tarn')]);
colormap(cmocean('tarn'));
% clim([0 4])
set(gca, 'ZDir','reverse')
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
annotation('line','LineStyle','--','Color',black,'LineWidth',0.75,'Position', [0.06 0.56 0.165 0.32]);% %%% Observed velocity
annotation('line','LineStyle','--','Color',black,'LineWidth',0.75,'Position', [0.31 0.75 -0.08 0.161]);% %%% Observed velocity
annotation('ellipse',[0.225 0.875 0.021 0.01],'Color',orange,'LineWidth',2,'rotation',70);
% ylim([min(ele_lat) 59.2])
freezeColors;



%%% Observed velocity
ax3 = subplot('position',[0.38 0.555 0.25 0.4]);
annotation('textbox',[0.375 0.99 0.15 0.01],'String','C','FontSize',fontsize+2,'fontweight','normal','LineStyle','None');

%%% Linear-fit velocity
ax4 = subplot('position',[0.73 0.555 0.25 0.4]);
annotation('textbox',[0.725 0.99 0.15 0.01],'String','D','FontSize',fontsize+2,'fontweight','normal','LineStyle','None');


%------ Data for plotting temperature
addpath ../observations/
temp = ncread('mavs2_20210718_level1.nc','__xarray_dataarray_variable__');
depth_temp = ncread('mavs2_20210718_level1.nc','depth')/1000;
time_temp = double(ncread('mavs2_20210718_level1.nc','time'));
time_temp = time_temp'/3600;%%% in hours
time_temp(end)=48;
% time_temp = time_temp/24;%%% in days
plot_tidx = 1:20:length(time_temp);
meanT = mean(temp,'all','omitnan');


%%% Temperature
ax5 = subplot('position',[0.06 0.12 0.26 0.33]);
annotation('textbox',[0.025 0.475 0.15 0.01],'String','E','FontSize',fontsize+2,'fontweight','normal','LineStyle','None');
pcolor(time_temp(plot_tidx),depth_temp,temp(plot_tidx,:)');shading flat;
hold on;
contour(time_temp(plot_tidx),depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
hold off;
xlabel('Time (hours)')
ylabel('Depth (km)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([meanT-1 meanT+2])
colormap(cmocean('balance'))
title('Conservative temperature','Fontsize',fontsize+2);
xticks([0:6:48])
h5=colorbar(ax5);
set(h5,'Position',[0.325 0.125 0.007 0.3]);
set(get(h5,'Title'),'String','   (^oC)');


%%% Reconstructed dbdz using the observed velocity
ax6 = subplot('position',[0.38 0.055 0.25 0.4]);
annotation('textbox',[0.375 0.49 0.15 0.01],'String','F','FontSize',fontsize+2,'fontweight','normal','LineStyle','None');

%%% Reconstructed dbdz using the linear-fit velocity
ax7 = subplot('position',[0.73 0.055 0.25 0.4]);
annotation('textbox',[0.725 0.49 0.15 0.01],'String','G','FontSize',fontsize+2,'fontweight','normal','LineStyle','None');
