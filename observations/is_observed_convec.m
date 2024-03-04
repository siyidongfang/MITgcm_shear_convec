
clear;close all;
addpath ../instability/
load_colors;
showplot = false;
fontsize = 20;

%%% Velocity data:
load('MAVS2_velocity.mat','uu','vv','ww','depth','time','topo');

time_adcp = double(ncread('MAVS2_24606.nc','time')) + 33.000323*1e6; % microseconds since 2021-07-07 06:00:33.000323
Time_adcp = datetime(2021,7,7) + hours(round(time_adcp/1e6/3600)+6);
Time_adcp.Format = 'yyyy-MM-dd';

% Find time indices in the ADCP data of two tidal cycles:
nStart = 1031;
nEnd = 1222;
tidx = nStart:nEnd;
zidx = 4:18;
uselect = uu(zidx,tidx)';
vselect = vv(zidx,tidx)';
wselect = ww(zidx,tidx)';
time_uw = time(tidx)/3600; %%% in hours
time_uw = time_uw-time_uw(1);

depth_uw = depth(zidx);

%%% Temperature data of two tidal cycles:
temp = ncread('mavs2_20210718_level1.nc','__xarray_dataarray_variable__');
depth_temp = ncread('mavs2_20210718_level1.nc','depth');
time_temp = double(ncread('mavs2_20210718_level1.nc','time'));

meanT = mean(temp,'all','omitnan');
T_tavg = mean(temp,'omitnan');
tAlpha = 2e-4;
gravity = 9.81;

dT_dz = diff(temp,1,2)./repmat(diff(-depth_temp'),[length(temp) 1]);
N2 = gravity*tAlpha*dT_dz;
dTavg_dz = diff(T_tavg)./diff(-depth_temp');

N2avg = gravity*tAlpha*dTavg_dz;
meanN2 = mean(N2avg);

depth_n2 = 0.5*(depth_temp(1:end-1)+depth_temp(2:end))';
time_temp = time_temp'/3600;%%% in hours


%%% Interpolate velocity and stratification on the same 
%%% vertical grid and time grid
depth_uw2 = depth(min(zidx)-1:max(zidx)+1);
depth_uw_stagger = 0.5*(depth_uw2(1:end-1)+depth_uw2(2:end));

[DD_temp,TT_temp] = meshgrid(depth_temp,time_temp);
[DD_uw_stagger,TT_uw_stagger] = meshgrid(depth_uw_stagger,time_uw);

temp_uw_stagger = interp2(DD_temp,TT_temp,temp,DD_uw_stagger,TT_uw_stagger);
N2_uwgrid = gravity*tAlpha.*diff(temp_uw_stagger,1,2)./repmat(diff(-depth_uw_stagger),[length(temp_uw_stagger) 1]);

figure(6)
clf;set(gcf,'Color','w');
subplot(1,2,1)
pcolor(time_uw,depth_uw_stagger,temp_uw_stagger');shading flat;colorbar;
hold on;
contour(time_uw,depth_uw_stagger,temp_uw_stagger',meanT-2:0.5:meanT+2,'Color',black);
hold off;
xlabel('Time (hours)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([meanT-2 meanT+2])
title('Temperature (^oC)')

subplot(1,2,2)
pcolor(time_uw,depth_uw,N2_uwgrid');shading flat;colorbar;
hold on;
contour(time_uw,depth_uw_stagger,temp_uw_stagger',meanT-2:0.5:meanT+2,'Color',black);
hold off;
xlabel('Time (hours)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
title('N^2 (1/s^2)')
clim([-3 3]/1e5)
colormap(cmocean('balance'));


%% time-mean, depth-averaged N^2 + observed velocities
adv1 = uselect*meanN2*sind(topo)+wselect*meanN2*cosd(topo);


%% time-mean, depth-varying N^2 + observed velocities
adv2 = uselect.*N2avg*sind(topo)+wselect.*N2avg*cosd(topo);


%% observed N^2 + observed velocities
adv3 = uselect.*N2*sind(topo)+wselect.*N2*cosd(topo);


if(showplot)
figure(1);
clf;set(gcf,'Color','w');
subplot(1,2,1)
pcolor(time_temp,depth_temp,temp');shading flat;colorbar;
hold on;
contour(time_temp,depth_temp,temp',[meanT-2:0.5:meanT+2],'Color',black);
hold off;
xlabel('Time (hours)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([meanT-2 meanT+2])
title('Temperature (^oC)')

subplot(1,2,2)
pcolor(time_temp,depth_n2,N2');shading flat;colorbar;
hold on;
contour(time_temp,depth_temp,temp',meanT-2:0.5:meanT+2,'Color',black);
hold off;
xlabel('Time (hours)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
title('N^2 (1/s^2)')
clim([-3 3]/1e5)
colormap(cmocean('balance'));

figure(2);
clf;set(gcf,'Color','w');
subplot(1,2,1)
plot(T_tavg,depth_temp,'LineWidth',2);
axis ij;set(gca,'Fontsize',fontsize);
grid on;grid minor;
ylabel('Depth (m)')
title('2-day mean temperature (^oC)')

subplot(1,2,2)
plot(N2avg,depth_n2,'LineWidth',2);
axis ij;set(gca,'Fontsize',fontsize);
grid on;grid minor;
ylabel('Depth (m)')
title('2-day mean N^2 (1/s^2)')

figure(4);
clf;set(gcf,'Color','w');
set(gcf,'Position', [93 321 1680 364])
subplot(1,3,1)
pcolor(time_uw,depth_uw,uselect');shading interp;colorbar;colormap(redblue)
xlabel('Time (hours)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([-0.3 0.3])
title('u (m/s)')

subplot(1,3,2)
pcolor(time_uw,depth_uw,vselect');shading interp;colorbar;colormap(redblue)
xlabel('Time (hours)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([-0.3 0.3])
title('v (m/s)')

subplot(1,3,3)
pcolor(time_uw,depth_uw,wselect');shading interp;colorbar;colormap(redblue)
xlabel('Time (hours)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([-0.3 0.3]/6)
title('w (m/s)')
end



