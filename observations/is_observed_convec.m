
clear;close all;
addpath ../instability/

temp = ncread('mavs2_20210718_level1.nc','__xarray_dataarray_variable__');
depth = ncread('mavs2_20210718_level1.nc','depth');
time = double(ncread('mavs2_20210718_level1.nc','time'));

time_adcp = double(ncread('MAVS2_24606.nc','time')) + 33.000323*1e6; % microseconds since 2021-07-07 06:00:33.000323
Time_adcp = datetime(2021,7,7) + hours(round(time_adcp/1e6/3600)+6);
Time_adcp.Format = 'yyyy-MM-dd';


u = ncread('MAVS2_24606.nc','u');
v = ncread('MAVS2_24606.nc','v');
w = ncread('MAVS2_24606.nc','w');
depth_adcp = ncread('MAVS2_24606.nc','depth');


%%% Select velocity, temperature of two tidal cycles
%%
fontsize = 20;

meanT = mean(temp,'all','omitnan');
T_tavg = mean(temp,'omitnan');
tAlpha = 2e-4;
gravity = 9.81;

dTavg_dz = (T_tavg(1:end-1)-T_tavg(2:end))./(-(depth(1:end-1)-depth(2:end)))';

diff(T_tavg)./diff(depth');
N2avg = gravity*tAlpha*dTavg_dz;
meanN = mean(N2avg);

figure(2)
plot(T_tavg,depth);axis ij;
figure(3)
plot(N2avg,0.5*(depth(1:end-1)+depth(2:end)));axis ij;


figure(1);clf;set(gcf,'Color','w');
pcolor(time/3600,depth,temp');shading flat;colorbar;colormap(redblue);
xlabel('Time (hours)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([meanT-2 meanT+2])

%%% Find time indices in the ADCP data
nStart = 1031;
nEnd = 1222;
tidx = nStart:nEnd;
uu = u(tidx,:);
vv = v(tidx,:);
ww = w(tidx,:);
tt = time_adcp(tidx)/1e6/3600; %%% in hours
tt = tt-tt(1)+15/60;


%%
figure(4);
clf;set(gcf,'Color','w');
set(gcf,'Position', [93 321 1680 364])
subplot(1,3,1)
pcolor(tt,depth_adcp,uu');shading interp;colorbar;colormap(redblue)
xlabel('Time (hours)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([-0.3 0.3])
title('u (m/s)')

subplot(1,3,2)
pcolor(tt,depth_adcp,vv');shading interp;colorbar;colormap(redblue)
xlabel('Time (hours)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([-0.3 0.3])
title('v (m/s)')


subplot(1,3,3)
pcolor(tt,depth_adcp,ww');shading interp;colorbar;colormap(redblue)
xlabel('Time (hours)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([-0.3 0.3]/6)
title('w (m/s)')





