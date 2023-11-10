%%%
%%% velocitySpectra.m
%%%
%%% Plot the spectra of velocity and shear for mooring data

addpath /Users/ysi/MITgcm_shear_convec/analysis/colormaps/
addpath /Users/ysi/MITgcm_shear_convec/analysis/colormaps/cmocean/
addpath /Users/ysi/MITgcm_shear_convec/observations/moorings_gridded_adcp_v2
ncfname = 'MAVS2_24606.nc';

tt = ncread(ncfname,'temperature')';
u_meridional = ncread(ncfname,'u')';
v_zonal = ncread(ncfname,'v')';
ww = ncread(ncfname,'w')'; 
pp = ncread(ncfname,'pressure')'; 
depth = ncread(ncfname,'depth')'; 
time  = double(ncread(ncfname,'time'))'/1e6; %%% convert microsecond to second

%%% Calculate along-canyon velocity uu and cross-slope velocity vv
dy_canyon = 
dx_canyon = 
angle_canyon = atand(dy_canyon/dx_canyon); %%% The angle between the zonal direction and the canyon
uu = u_meridional*cosd(angle_canyon)-v_zonal*sind(angle_canyon);
vv = u_meridional*sind(angle_canyon)+v_zonal*cosd(angle_canyon);



% dz = diff(depth); %%% DOUBLE CHECK dz!!
% dz = dz(1)
% dudz = diff(uu)./dz;
% dvdz = diff(vv)./dz;
% dwdz = diff(ww)./dz;


figure(1)
subplot(2,1,1)
plot(time/86400,pp);ylim([1150 1170]);grid on;grid minor;
title('Pressure (dbar)');xlabel('time (days)')
subplot(2,1,2)
plot(time/86400,tt);grid on;grid minor;
title('Temperature (^oC)');xlabel('time (days)')

t_idx = 97:96*10;

figure(2)
subplot(3,1,1)
pcolor(u_meridional(:,t_idx));shading interp;colorbar;colormap(redblue);clim([-0.3 0.3])
subplot(3,1,2)
pcolor(v_zonal(:,t_idx));shading interp;colorbar;colormap(redblue);clim([-0.3 0.3])
subplot(3,1,3)
pcolor(ww(:,t_idx));shading interp;colorbar;colormap(redblue);clim([-0.3 0.3]/10)






