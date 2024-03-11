

clear;close all;

addpath /Users/ysi/MITgcm_shear_convec/observations/topography/
addpath /Users/ysi/MITgcm_shear_convec/analysis/colormaps/
addpath /Users/ysi/MITgcm_shear_convec/analysis/colormaps/cmocean/
addpath /Users/ysi/MITgcm_shear_convec/observations/moorings_gridded_adcp_v2/

%%% Calculate along-canyon velocity uu and cross-slope velocity vv
load("topog1D_28km.mat",'lat9','lon9','depth9')

ncfname = 'MAVS2_24606.nc';
lat = 54.18229903565569; 
lat1 = lat9(9);
lon1 = lon9(9);
lat2 = lat9(10);
lon2 = lon9(10);
t_idx = 32:8574;
depth1 = depth9(9);
depth2 = depth9(10);
zidx = 4:18;


% ncfname = 'MAVS1_24608.nc';
% lat = 54.197478060955085; 
% lat1 = lat9(8);
% lon1 = lon9(8);
% lat2 = lat9(9);
% lon2 = lon9(9);
% t_idx = 40:8756;
% depth1 = depth9(8);
% depth2 = depth9(9);
% zidx = 5:18;


% ncfname = 'MP1_24839.nc';
% lat = 54.23853400073004;
% lat1 = lat9(5);
% lon1 = lon9(5);
% lat2 = lat9(6);
% lon2 = lon9(6);
% t_idx = 12:644;
% depth1 = depth9(5);
% depth2 = depth9(6);
% zidx = 3:36;

distance12 = distance(lat2,lon2,lat1,lon1,referenceEllipsoid('GRS80','m'));
topo = atand(-(depth2-depth1)./distance12);

tt = ncread(ncfname,'temperature')';
u_meridional = ncread(ncfname,'u')';
v_zonal = ncread(ncfname,'v')';
pp = ncread(ncfname,'pressure')'; 
depth = ncread(ncfname,'depth')'; 
time  = double(ncread(ncfname,'time'))'/1e6; %%% convert microsecond to second

dy_canyon = distance(lat1,lon2,lat2,lon2,referenceEllipsoid('GRS80','km'));
dx_canyon = distance(lat2,lon1,lat2,lon2,referenceEllipsoid('GRS80','km'));

angle_canyon = atand(dy_canyon/dx_canyon); %%% The angle between the zonal direction and the canyon

%%% Velocity in the normal coordinate
uu_tilde = u_meridional*cosd(angle_canyon)-v_zonal*sind(angle_canyon);
vv_tilde = u_meridional*sind(angle_canyon)+v_zonal*cosd(angle_canyon);
ww_tilde = ncread(ncfname,'w')'; 

%%% Velocity in the slope-aligned coordinate
uu = uu_tilde*cosd(topo)+ww_tilde*sind(topo);
ww = -uu_tilde*sind(topo)+ww_tilde*cosd(topo);
vv = vv_tilde;


fontsize = 20;
% t_idx_plot = 97:96*10;
t_idx_plot = 1:length(time);

dz = diff(depth); %%% DOUBLE CHECK dz!!
dz = dz(1)
dudz = diff(uu)./dz;
% dvdz = diff(vv)./dz;
% dwdz = diff(ww)./dz;

save('MAVS2_velocity.mat',...
    'uu','vv','ww',...
    'uu_tilde','vv_tilde','ww_tilde',...
    'u_meridional','v_zonal',...
    'depth','angle_canyon','time','topo')

% figure(1);clf;set(gcf,'color','w');
% subplot(2,1,1)
% plot(time/86400,pp);ylim([1150 1170]);grid on;grid minor;
% title('Pressure (dbar)');xlabel('time (days)')
% set(gca,'FontSize',fontsize)
% subplot(2,1,2)
% plot(time/86400,tt);grid on;grid minor;
% title('Temperature (^oC)');xlabel('time (days)')
% set(gca,'FontSize',fontsize)

figure(1);
clf;set(gcf,'color','w','Position',[249 367 931 599]);
subplot(3,1,1)
pcolor(time(t_idx_plot)/86400,depth(zidx),u_meridional(zidx,t_idx_plot));
shading interp;colorbar;colormap(redblue);clim([-0.4 0.4])
title('Zonal velocity (m/s)')
ylabel('Depth (m)');xlabel('time (days)');set(gca,'FontSize',fontsize)
axis ij;
subplot(3,1,2)
pcolor(time(t_idx_plot)/86400,depth(zidx),v_zonal(zidx,t_idx_plot));
shading interp;colorbar;colormap(redblue);clim([-0.4 0.4])
title('Meridional velocity (m/s)')
ylabel('Depth (m)');xlabel('time (days)');set(gca,'FontSize',fontsize)
axis ij;
subplot(3,1,3)
pcolor(time(t_idx_plot)/86400,depth(zidx),ww_tilde(zidx,t_idx_plot));
shading interp;colorbar;colormap(redblue);clim([-0.05 0.05])
title('Vertical velocity (m/s)')
ylabel('Depth (m)');xlabel('time (days)');set(gca,'FontSize',fontsize)
axis ij;


figure(2);
clf;set(gcf,'color','w','Position',[249 367 931 599]);
subplot(3,1,1)
pcolor(time(t_idx_plot)/86400,depth(zidx),uu_tilde(zidx,t_idx_plot));
shading interp;colorbar;colormap(redblue);clim([-0.4 0.4])
title('Up-canyon velocity, normal coordinate (m/s)')
ylabel('Depth (m)');xlabel('time (days)');set(gca,'FontSize',fontsize)
axis ij;
subplot(3,1,2)
pcolor(time(t_idx_plot)/86400,depth(zidx),vv_tilde(zidx,t_idx_plot));
shading interp;colorbar;colormap(redblue);clim([-0.4 0.4]/2)
title('Cross-canyon velocity (m/s)')
ylabel('Depth (m)');xlabel('time (days)');set(gca,'FontSize',fontsize)
axis ij;
subplot(3,1,3)
pcolor(time(t_idx_plot)/86400,depth(zidx),ww_tilde(zidx,t_idx_plot));
shading interp;colorbar;colormap(redblue);clim([-0.05 0.05])
title('Vertical velocity, normal coordinate (m/s)')
ylabel('Depth (m)');xlabel('time (days)');set(gca,'FontSize',fontsize)
axis ij;



figure(3);
clf;set(gcf,'color','w','Position',[249 367 931 599]);
subplot(3,1,1)
pcolor(time(t_idx_plot)/86400,depth(zidx),uu(zidx,t_idx_plot));
shading interp;colorbar;colormap(redblue);clim([-0.4 0.4])
title('Up-canyon velocity, slope-alighed coordinate (m/s)')
ylabel('Depth (m)');xlabel('time (days)');set(gca,'FontSize',fontsize)
axis ij;
subplot(3,1,2)
pcolor(time(t_idx_plot)/86400,depth(zidx),vv(zidx,t_idx_plot));
shading interp;colorbar;colormap(redblue);clim([-0.4 0.4]/2)
title('Cross-canyon velocity (m/s)')
ylabel('Depth (m)');xlabel('time (days)');set(gca,'FontSize',fontsize)
axis ij;
subplot(3,1,3)
pcolor(time(t_idx_plot)/86400,depth(zidx),ww(zidx,t_idx_plot));
shading interp;colorbar;colormap(redblue);clim([-0.05 0.05])
title('Vertical velocity, slope-alighed coordinate (m/s)')
ylabel('Depth (m)');xlabel('time (days)');set(gca,'FontSize',fontsize)
axis ij;

