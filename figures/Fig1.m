
clear;close all;
addpath ../analysis/colormaps/
fontsize = 17;
load_colors;


%%%------ Data for plotting bathymetry
addpath ../observations/topography/
ncfname = 'blt_canyon_mb_qc.nc';
lat = ncread(ncfname,'lat')';
lon = ncread(ncfname,'lon')';
z = ncread(ncfname,'z');
load Rockall_gebco.mat

load('fig1/fig1.mat')
plot_tidx = 1:1:length(time_temp);

zz_rec = -depth_n2 + depth_n2(end);
buoy_init = repmat(meanN2*cosd(topo)*zz_rec,[length(buoy1fit) 1]);
buoy_fit = buoy1fit + buoy_init;
buoy_obs = buoy1obs + buoy_init;
n2_1obs = meanN2 + cosd(topo)*(n2_1obs-meanN2*cosd(topo));
n2_1fit = meanN2 + cosd(topo)*(n2_1fit-meanN2*cosd(topo));


figure(1)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 1400 700]);

%%% Canyon bathymetry
ax2 = subplot('position',[0.04 0.53 0.27 0.25]);
annotation('textbox',[0.025 0.993 0.15 0.01],'String','a','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
surf(lon,lat,z'/1000,'EdgeColor','None');
shading flat;
set(gca, 'ZDir','reverse')
clim([min(min(z))/1000-max(max(z))/1000 max(max(z))/1000])
xlabel('Longitude','interpreter','latex','Position',[-12.1 54.065 4.2]);
ylabel('Latitude','interpreter','latex','Position',[-12.1 54.2286 3.55])
zlabel('Depth (km)','interpreter','latex');
set(gca,'FontSize',fontsize,'TickDir', 'in','TickLength',[0.1, 0.014]);
title('Canyon','FontSize',fontsize+5,'interpreter','latex','Position',[-12.24 54.22 -1.4]);
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
mycolor=flip(cmocean('diff'));
% mycolor=cmocean('topo','negative');
mycolor=mycolor([53:5:128 129:1:end],:);
colormap(mycolor);
% colormap(cmocean('topo','negative'));
% clim([-3.5 3.5])
clim([-0.5 4])
freezeColors;
%%
%%% Rockall Trough
ax1 = subplot('position',[0.028 0.8 0.28 0.17]);
annotation('textbox',[0.025 0.76 0.15 0.01],'String','b','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
surf(ele_lon,ele_lat,-elev'/1000,'EdgeColor','None');
hold on;
contour3(ele_lon,ele_lat,-elev'/1000,[0 0],'EdgeColor',black);
% contour3(ele_lon,ele_lat,-elev'/1000,[0.45 0.45],'-.','EdgeColor',black,'ShowText','off');
hold off;
set(gca, 'ZDir','reverse')
shading flat;
% colormap(flip(cmocean('topo')));
% colormap(cmocean('topo','negative'));
% clim([-3.5 3.5])
clim([-0.5 4])
% colormap(flip(haxby))
% colormap([[0.1 0.1 0.1]*2;cmocean('tarn')]);
colormap(mycolor);
% clim([0 4])
set(gca, 'ZDir','reverse')
% view([-78 58]);
view([-84.98   52.81]);
xlabel('Longitude','interpreter','latex','Position',[-18.3013 49.6 7]);
ylabel('Latitude','interpreter','latex','Position',[-22 54.7165 1.5]);
zlabel('Depth (km)','interpreter','latex');
set(get(gca,'xlabel'),'rotation',67);
set(get(gca,'ylabel'),'rotation',-1.5);
axis tight;box off;
set(gca,'FontSize',fontsize,'TickDir', 'in','TickLength',[0.1, 0.014]);
title('Rockall Trough','FontSize',fontsize+5,'interpreter','latex','Position',[-5.5 57.5 0]);
xticks([-19:3:-6]);
annotation('line','LineStyle','--','Color',black,'LineWidth',0.75,'Position', [0.04 0.543 0.179 0.332]);
annotation('line','LineStyle','--','Color',black,'LineWidth',0.75,'Position', [0.31 0.765 -0.084 0.152]);
annotation('ellipse',[0.22 0.875 0.021 0.006],'Color',yellow,'LineWidth',2,'rotation',70);
% ylim([min(ele_lat) 59.2])
freezeColors;


%% Temperature

ax3 = subplot('position',[0.048 0.07 0.25 0.38]);
annotation('textbox',[0.028 0.482 0.15 0.01],'String','c','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(time_temp(plot_tidx)*24,depth_temp,temp(plot_tidx,:)');shading flat;
hold on;
contour(time_temp(plot_tidx)*24,depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
contour(time_temp(plot_tidx)*24,depth_n2,smooth_N2(plot_tidx,:)',[0 0],'Color',cyan,'LineWidth',0.75);
hold off;
xlabel('Time (hours)','interpreter','latex');
ylabel('Depth (m)','interpreter','latex');
set(gca,'Fontsize',fontsize);
axis ij;
% clim([meanT-1 meanT+2]);
clim([3.85 6.5]);
xlim([0 48])
colormap(cmocean('balance'))
title('Potential temperature','Fontsize',fontsize+5,'interpreter','latex');
xticks([0:6:48])
h3=colorbar(ax3);
set(h3,'Position',[0.303 0.135 0.007 0.28]);
set(get(h3,'Title'),'String','$\ \ \ \ (^\circ \mathrm{C})$','Fontsize',fontsize,'interpreter','latex');
%%

%%% Observed velocity
ax4 = subplot('position',[0.38 0.58 0.25 0.38]);
annotation('textbox',[0.36 0.993 0.15 0.01],'String','d','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(time_u*24,depth_u,uobs');
hold on;
contour(time_temp(plot_tidx)*24,depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
hold off;
shading flat;
xlabel('Time (hours)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (m)','interpreter','latex','FontSize',fontsize);
set(gca,'Fontsize',fontsize);
axis ij;
clim([-0.4 0.4])
colormap(cmocean('balance'))
title('Observed velocity $u_\mathrm{obs}$','Fontsize',fontsize+5,'interpreter','latex');
xticks([0:6:48])
h4=colorbar(ax4);
ylim([min(depth_temp) max(depth_temp)]);xlim([0 48])
set(h4,'Position',[0.635  0.645 0.007 0.28]);
set(get(h4,'Title'),'String','$\ \ \ \ (\mathrm{m/s})$','Fontsize',fontsize,'interpreter','latex');

%%% Linear-fit velocity
ax5 = subplot('position',[0.709 0.58 0.25 0.38]);
annotation('textbox',[0.69 0.993 0.15 0.01],'String','e','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(time_u*24,depth_u,ufit');
hold on;
contour(time_temp(plot_tidx)*24,depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
hold off;
shading flat;
xlabel('Time (hours)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (m)','interpreter','latex','FontSize',fontsize);
set(gca,'Fontsize',fontsize);
axis ij;
clim([-0.4 0.4])
colormap(cmocean('balance'))
title('Linear-fit velocity $U_\mathrm{fit}$','Fontsize',fontsize+5,'interpreter','latex');
xticks([0:6:48])
ylim([min(depth_temp) max(depth_temp)]);xlim([0 48])
h5=colorbar(ax5);
set(h5,'Position',[0.964  0.645 0.007 0.28]);
set(get(h5,'Title'),'String','$\ \ \ \ (\mathrm{m/s})$','Fontsize',fontsize,'interpreter','latex');


%%% Reconstructed dbdz using the observed velocity
ax6 = subplot('position',[0.38 0.07 0.25 0.38]);
annotation('textbox',[0.36 0.482 0.15 0.01],'String','f','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(time_temp(plot_tidx)*24,depth_n2,buoy_obs(plot_tidx,:)');
% pcolor(time_temp(plot_tidx)*24,depth_reconst_n,n2_1obs(plot_tidx,:)');
shading flat;
hold on;
contour(time_temp(plot_tidx)*24,depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
contour(time_temp(plot_tidx)*24,depth_reconst_n,n2_1obs(plot_tidx,:)',[0 0],'Color','c','LineWidth',2);
% contour(time_temp(plot_tidx)*24,depth_reconst_n,n2_1obs(plot_tidx,:)',[1 2 3 4]*1e-6,'Color','w','LineWidth',1.5);
% contour(time_temp(plot_tidx)*24,depth_n2,buoy_obs(plot_tidx,:)','Color',black);
hold off;
xlabel('Time (hours)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (m)','interpreter','latex','FontSize',fontsize);
set(gca,'Fontsize',fontsize);
axis ij;
% clim([-1.75 1.75]/1e5)
clim([-0.2 1.7]/1e3)
colormap(cmocean('balance'))
title('Reconstruct $\mathcal B_\mathrm{rec}$ using $u_\mathrm{obs}$','Fontsize',fontsize+5,'interpreter','latex');
% title('Reconstruct $\partial_{\tilde z} \mathcal B$ using $u_\mathrm{obs}$','Fontsize',fontsize+5,'interpreter','latex');
xticks([0:6:48])
h6=colorbar(ax6);
set(h6,'Position',[0.635 0.135 0.007 0.28]);
set(get(h6,'Title'),'String',{'$\ \ \ \ (\mathrm{m/s^2})$',''},'Fontsize',fontsize,'interpreter','latex');
% set(get(h6,'Title'),'String',{'$\ \ \ \ (1/\mathrm{s}^2)$',''},'Fontsize',fontsize,'interpreter','latex');
xlim([0 48])


%%% Reconstructed dbdz using the linear-fit velocity
ax7 = subplot('position',[0.709 0.07 0.25 0.38]);
annotation('textbox',[0.69 0.482 0.15 0.01],'String','g','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(time_temp(plot_tidx)*24,depth_n2,buoy_fit(plot_tidx,:)');
shading flat;
% pcolor(time_temp(plot_tidx)*24,depth_reconst_n,n2_1fit(plot_tidx,:)');shading interp;
hold on;
contour(time_temp(plot_tidx)*24,depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
contour(time_temp(plot_tidx)*24,depth_reconst_n,n2_1fit(plot_tidx,:)',[0 0],'Color','c','LineWidth',2);
% contour(time_temp(plot_tidx)*24,depth_reconst_n,n2_1fit(plot_tidx,:)',[1 2 3 4]*1e-6,'Color','w','LineWidth',1.5);
% contour(time_temp(plot_tidx)*24,depth_n2,buoy_fit(plot_tidx,:)','Color',black);
hold off;
xlabel('Time (hours)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (m)','interpreter','latex','FontSize',fontsize);
set(gca,'Fontsize',fontsize);
axis ij;
% clim([-1.75 1.75]/1e5)
clim([-0.2 1.7]/1e3)
colormap(cmocean('balance'))
title('Reconstruct $\mathcal B_\mathrm{rec}$ using $U_\mathrm{fit}$','Fontsize',fontsize+5,'interpreter','latex');
% title('Reconstruct $\partial_{\tilde z} \mathcal B$ using $U_\mathrm{fit}$','Fontsize',fontsize+5,'interpreter','latex');
xticks([0:6:48])
h7=colorbar(ax7);
set(h7,'Position',[0.964 0.135 0.007 0.28]);
set(get(h7,'Title'),'String',{'$\ \ \ \ (\mathrm{m/s^2})$',''},'Fontsize',fontsize,'interpreter','latex');
% set(get(h7,'Title'),'String',{'$\ \ \ \ (1/\mathrm{s}^2)$',''},'Fontsize',fontsize,'interpreter','latex');
xlim([0 48])

print('-dpng','-r300','fig1/fig1_matlab.png');
