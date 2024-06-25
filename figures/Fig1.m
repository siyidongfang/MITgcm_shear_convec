
clear;close all;
addpath ../analysis/colormaps/
addpath freezeColors/
fontsize = 16;
load_colors;

calc_Fig1;

figure(1)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 1400 700]);

%%% Canyon bathymetry
ax2 = subplot('position',[0.04 0.53 0.27 0.25]);
annotation('textbox',[0.025 0.99 0.15 0.01],'String','A','FontSize',fontsize+3,'fontweight','normal','LineStyle','None');
surf(lon,lat,z'/1000,'EdgeColor','None');
shading flat;
set(gca, 'ZDir','reverse')
clim([min(min(z))/1000-max(max(z))/1000 max(max(z))/1000])
set(gca,'FontSize',fontsize,'TickDir', 'in','TickLength',[0.1, 0.014]);
xlabel('Longitude','Position',[-12.1 54.065 4]);
ylabel('Latitude','Position',[-12.1 54.2286 3.45])
zlabel('Depth (km)');
title('Canyon','FontSize',fontsize+3,'Position',[-12.24 54.22 -1.4]);
% view([-78 58]);
view([-84.98 52.81]);
set(get(gca,'xlabel'),'rotation',77);
set(get(gca,'ylabel'),'rotation',-1);
axis tight;
ax = gca;
ax.GridColor = [0.1, 0.1, 0.1]*5; 
box off
xticks([-12.2:0.1:-11.8]);
% ylim([54.1 54.38])
% colormap([[0.1 0.1 0.1]*2;cmocean('rain')]);
mycolor=cmocean('tarn');
% mycolor=mycolor([80:3:128 129:1:end],:);
colormap(mycolor);
% clim([-0.5 4])
% colormap(cmocean('topo','negative'));
clim([-3.5 3.5])
freezeColors;

%%% Rockall Trough
ax1 = subplot('position',[0.028 0.8 0.28 0.17]);
annotation('textbox',[0.025 0.74 0.15 0.01],'String','B','FontSize',fontsize+3,'fontweight','normal','LineStyle','None');
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
% clim([-0.5 4])
% colormap(flip(haxby))
% colormap([[0.1 0.1 0.1]*2;cmocean('tarn')]);
colormap(mycolor);
% clim([0 4])
set(gca, 'ZDir','reverse')
% view([-78 58]);
view([-84.98   52.81]);
zlabel('(km)')
xlabel('Longitude','Position',[-18.3013 49.6 6]);
ylabel('Latitude','Position',[-22 54.7165 1]);
zlabel('Depth (km)')
set(get(gca,'xlabel'),'rotation',67);
set(get(gca,'ylabel'),'rotation',-1.5);
axis tight;box off;
set(gca,'FontSize',fontsize,'TickDir', 'in','TickLength',[0.1, 0.014]);
title('Rockall Trough','FontSize',fontsize+3,'Position',[-5.5 57.5 -0.4]);
xticks([-19:3:-6]);
annotation('line','LineStyle','--','Color',black,'LineWidth',0.75,'Position', [0.04 0.543 0.178 0.337]);
annotation('line','LineStyle','--','Color',black,'LineWidth',0.75,'Position', [0.31 0.765 -0.085 0.153]);
annotation('ellipse',[0.22 0.875 0.021 0.01],'Color',orange,'LineWidth',2,'rotation',70);
% ylim([min(ele_lat) 59.2])
freezeColors;


%%% Temperature
ax3 = subplot('position',[0.045 0.07 0.25 0.38]);
annotation('textbox',[0.025 0.482 0.15 0.01],'String','C','FontSize',fontsize+3,'fontweight','normal','LineStyle','None');
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
title('Conservative temperature','Fontsize',fontsize+3);
xticks([0:6:48])
h3=colorbar(ax3);
set(h3,'Position',[0.3 0.12 0.007 0.3]);
set(get(h3,'Title'),'String','   (^oC)');

%%


%%% Observed velocity
ax4 = subplot('position',[0.38 0.58 0.25 0.38]);
annotation('textbox',[0.375 0.99 0.15 0.01],'String','D','FontSize',fontsize+3,'fontweight','normal','LineStyle','None');
pcolor(time_temp(plot_tidx),depth_temp,temp(plot_tidx,:)');shading flat;
xlabel('Time (hours)')
ylabel('Depth (km)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([meanT-1 meanT+2])
colormap(cmocean('balance'))
title('Observed velocity','Fontsize',fontsize+3);
xticks([0:6:48])
h4=colorbar(ax4);
set(h4,'Position',[0.635  0.63 0.007 0.3]);
set(get(h4,'Title'),'String','   (m/s)');


%%% Linear-fit velocity
ax5 = subplot('position',[0.715 0.58 0.25 0.38]);
annotation('textbox',[0.71 0.99 0.15 0.01],'String','E','FontSize',fontsize+3,'fontweight','normal','LineStyle','None');
pcolor(time_temp(plot_tidx),depth_temp,temp(plot_tidx,:)');shading flat;
xlabel('Time (hours)')
ylabel('Depth (km)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([meanT-1 meanT+2])
colormap(cmocean('balance'))
title('Linear-fit velocity','Fontsize',fontsize+3);
xticks([0:6:48])
h5=colorbar(ax5);
set(h5,'Position',[0.97  0.63 0.007 0.3]);
set(get(h5,'Title'),'String','   (m/s)');



%%% Reconstructed dbdz using the observed velocity
ax6 = subplot('position',[0.38 0.07 0.25 0.38]);
annotation('textbox',[0.375 0.482 0.15 0.01],'String','F','FontSize',fontsize+3,'fontweight','normal','LineStyle','None');
pcolor(time_temp(plot_tidx),depth_temp,temp(plot_tidx,:)');shading flat;
xlabel('Time (hours)')
ylabel('Depth (km)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([meanT-1 meanT+2])
colormap(cmocean('balance'))
title('Reconstructed \partial b/\partial z (observed)','Fontsize',fontsize+3);
xticks([0:6:48])
h6=colorbar(ax6);
set(h6,'Position',[0.635 0.12 0.007 0.3]);
set(get(h6,'Title'),'String','   (^oC)');


%%% Reconstructed dbdz using the linear-fit velocity
ax7 = subplot('position',[0.715 0.07 0.25 0.38]);
annotation('textbox',[0.71 0.482 0.15 0.01],'String','G','FontSize',fontsize+3,'fontweight','normal','LineStyle','None');
pcolor(time_temp(plot_tidx),depth_temp,temp(plot_tidx,:)');shading flat;
xlabel('Time (hours)')
ylabel('Depth (km)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([meanT-1 meanT+2])
colormap(cmocean('balance'))
title('Reconstructed \partial b/\partial z (linear-fit)','Fontsize',fontsize+3);
xticks([0:6:48])
h7=colorbar(ax7);
set(h7,'Position',[0.97 0.12 0.007 0.3]);
set(get(h7,'Title'),'String','   (m/s)');