%%%
%%% plot_RockallTopog.m
%%%
%%% Plot the Rockall Trough topography, and calculate and plot the topographic slope
%%% Determine whether the slope is critical slope for local tides

addpath /Users/ysi/MITgcm_BLT/analysis/observations/topography
addpath /Users/ysi/MITgcm_BLT/analysis/functions/

load Rockall_gebco.mat
load_colors;
elev(elev>=0)=NaN;

%%% Find the canyon
% lat_max = 54.25;
% lat_min = 54.165;
% lon_max = -11.81;
% lon_min = -11.96; 


lat_max = 54.25+0.25;
lat_min = 54.165-0.25;
lon_max = -11.81+0.25;
lon_min = -11.96-0.25; 

lat_idx = find((ele_lat<=lat_max).*(ele_lat>=lat_min)==1);
lon_idx = find((ele_lon<=lon_max).*(ele_lon>=lon_min)==1);
lat_canyon = ele_lat(lat_idx);
lon_canyon = ele_lon(lon_idx);
elev_canyon = elev(lon_idx,lat_idx);

%%% Make figures
fontsize = 16;

figure(10);
clf;set(gcf,'color','w');
subplot(1,2,1)
pcolor(ele_lon,ele_lat,elev')
shading flat;colorbar;colormap(cmocean('rain'));
set(gca,'FontSize',fontsize);
title('Rockall Trough topography (m)');
xlabel('Longitude'); ylabel('Latitude')
clim([-4500 -300])

subplot(1,2,2)
pcolor(lon_canyon,lat_canyon,elev_canyon')
hold on;
contour(lon_canyon,lat_canyon,elev_canyon',[-800:50:-200])
shading flat;colorbar;
% colormap(cmocean('rain'));
set(gca,'FontSize',fontsize);
title('Canyon topography (m)');
xlabel('Longitude'); ylabel('Latitude')
% clim([-4500 -300])
%%% Select several cross-sections, and calculate the topographic slopes


figure(9)
contour(lon_canyon,lat_canyon,elev_canyon',[-2500:25:-200])
shading flat;colorbar;
shading flat;colorbar;
% colormap(cmocean('rain'));
set(gca,'FontSize',fontsize);
title('Canyon topography (m)');
xlabel('Longitude'); ylabel('Latitude')
grid on;grid minor;

 