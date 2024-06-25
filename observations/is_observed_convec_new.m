
clear;close all;
addpath ../instability/
addpath ../analysis/colormaps/
addpath /Users/ysi/Software/gsw_matlab_v3_06_11/
addpath /Users/ysi/Software/gsw_matlab_v3_06_11/library/
addpath /Users/ysi/Software/gsw_matlab_v3_06_11/thermodynamics_from_t/

load_colors;
showplot = true;
fontsize = 20;

%%% Velocity data:
load('MAVS2_LinearShear.mat');
load('MAVS2_velocity.mat')

time_adcp = double(ncread('MAVS2_24606.nc','time')) + 33.000323*1e6; % microseconds since 2021-07-07 06:00:33.000323
Time_adcp = datetime(2021,7,7) + hours(round(time_adcp/1e6/3600)+6);
Time_adcp.Format = 'yyyy-MM-dd';

% Find time indices in the ADCP data of two tidal cycles:
nStart = 1031-1;
nEnd = nStart+24*2*4+3;
% nStart = 263-1;
% nEnd = nStart+24*2*4+3;
tidx = nStart:nEnd;
zidx = 4:18;
ufit = u_reconstruct(zidx,tidx)';
uobs = uu(zidx,tidx)';

depth_u = depth(zidx);

meanUobs = mean(uobs,'omitnan');
meanUfit = mean(ufit,'omitnan');

time_u = time(tidx)/3600; %%% in hours
time_u = time_u-time_u(1);
time_u = time_u/24; %%% in days

Uzavgobs = mean((uobs),2);
Uzavgfit = mean((ufit),2);

figure(24)
clf;
plot(time_u,Uzavgobs/100);
hold on;
plot(time_u,shear_linear(tidx));
plot(time_u,Uzavgfit/100,'--');


figure(25)
plot(meanUobs,depth_u,'LineWidth',2);axis ij;
hold on;
plot(meanUfit,depth_u,'LineWidth',2);
grid on;
grid minor
set(gca,'Fontsize',fontsize);
set(gcf,'Color','w')
ylabel('Depth (m)')
title('Mean flow (m/s)')

uobs_detrend = uobs-mean(uobs);

ufit_detrend = ufit-mean(ufit);


%%

% Temperature data of two tidal cycles:
temp = ncread('mavs2_20210718_level1.nc','__xarray_dataarray_variable__');
depth_temp = ncread('mavs2_20210718_level1.nc','depth');
time_temp = double(ncread('mavs2_20210718_level1.nc','time'));

% addpath moorings_gridded_thermistors/level1/mavs2/
% temp = ncread('mavs2_20210710.nc','__xarray_dataarray_variable__');
% depth_temp = ncread('mavs2_20210710.nc','depth');
% time_temp = double(ncread('mavs2_20210710.nc','time'));

time_temp = time_temp'/3600;%%% in hours
time_temp = time_temp/24;%%% in days

meanT = mean(temp,'all','omitnan');
T_tavg = mean(temp,'omitnan');
dTavg_dz = diff(T_tavg)./diff(-depth_temp');

%%% Correct stratification using a linear relationship between T and S
o1 = 0.078611022276845;
o2 = 34.638914089176580;
salt = o1*temp+o2;
meanS = mean(salt,'all','omitnan');

%%% Use the GSW toolbox to calculate buoyancy frequency
lon_MAVS2 = -11.843;
lat_MAVS2 = 54.182;

[SA, in_ocean] = gsw_SA_from_SP(salt,depth_temp,lon_MAVS2,lat_MAVS2);
CT = gsw_CT_from_pt(SA,temp);

Nt_temp = length(time_temp);
Nz_temp = length(depth_temp);
N2 = NaN*zeros(length(time_temp),Nz_temp-1);
for ii = 1:Nt_temp
    [N2(ii,:), depth_n2] = gsw_Nsquared(SA(ii,:)',CT(ii,:)',depth_temp,lat_MAVS2);
end

depth_n2 = depth_n2';
N2_zavg = mean(N2,'omitnan');
meanN2 = mean(N2_zavg);

%%% Interpolate low-resolution velocity data onto the grid of high-res
%%% N2 grid
[DD_u,TT_u] = meshgrid(depth_u,time_u);
[DD_n2,TT_n2] = meshgrid(depth_n2,time_temp);

uobs_n2grid = interp2(DD_u,TT_u,uobs_detrend,DD_n2,TT_n2,'linear');
ufit_n2grid = interp2(DD_u,TT_u,ufit_detrend,DD_n2,TT_n2,'linear');

adv1obs = uobs_n2grid*meanN2*sind(topo);
adv1fit = ufit_n2grid*meanN2*sind(topo);

dt = 86400*(time_temp(2)-time_temp(1));

buoy1obs =  - cumsum(adv1obs*dt,1);
buoy1fit =  - cumsum(adv1fit*dt,1);

Nt = 2*(24*3600+50*60);

[n20tz, depth_reconst_n] = gsw_Nsquared(0.5*(SA(1,1:end-1)+SA(1,2:end)),...
    0.5*(CT(1,1:end-1)+CT(1,2:end)),depth_n2,lat_MAVS2);

n2_1obs = meanN2*cosd(topo) + diff(buoy1obs,1,2)./diff(-depth_n2);
n2_1fit = meanN2*cosd(topo) + diff(buoy1fit,1,2)./diff(-depth_n2);


%% Save data
save('fig1.mat','T_tavg','N2_zavg',...
'temp','N2','time_u','time_temp',...
'depth_temp','depth_n2','depth_u','depth_reconst_n',...
'uobs','ufit','n2_1obs','n2_1fit','meanT')



%% Make figures
addpath /Users/ysi/MITgcm_shear_convec/analysis/
load_colors;

plot_tidx = 1:10:length(time_temp);

figure(2);
clf;set(gcf,'Color','w','Position', [669 91 718 416]);
subplot(1,2,1)
plot(T_tavg,depth_temp,'LineWidth',2);
axis ij;set(gca,'Fontsize',fontsize);
grid on;grid minor;
ylabel('Depth (m)')
title('Time mean temperature (^oC)')

subplot(1,2,2)
plot(N2_zavg,depth_n2,'LineWidth',2);
axis ij;set(gca,'Fontsize',fontsize);
grid on;grid minor;
ylabel('Depth (m)')
title('Time mean N^2 (1/s^2)')

figure(4);
clf;set(gcf,'Color','w');
set(gcf,'Position', [56 352 1865 305]);
subplot(1,3,1)
pcolor(time_u,depth_u,uobs');shading interp;colorbar;colormap(redblue)
xlabel('Time (days)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([-0.3 0.3])
title('u (m/s)')

subplot(1,3,2)
pcolor(time_u,depth_u,uobs_detrend');shading flat;colorbar;colormap(redblue)
xlabel('Time (days)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([-0.3 0.3])
title('u (m/s): detrend')

subplot(1,3,3)
pcolor(time_temp,depth_n2,uobs_n2grid');shading interp;colorbar;colormap(redblue)
xlabel('Time (days)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([-0.3 0.3])
title('u (m/s): detrend, n2grid')



figure(5);
clf;set(gcf,'Color','w');
set(gcf,'Position', [56 352 1865 305]);
subplot(1,3,1)
pcolor(time_u,depth_u,ufit');shading flat;colorbar;colormap(redblue)
xlabel('Time (days)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([-0.3 0.3])
title('u (m/s): linear fit')

subplot(1,3,2)
pcolor(time_u,depth_u,ufit_detrend');shading flat;colorbar;colormap(redblue)
xlabel('Time (days)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([-0.3 0.3])
title('u (m/s): linear fit, detrend')


subplot(1,3,3)
pcolor(time_temp,depth_n2,ufit_n2grid');shading flat;colorbar;colormap(redblue)
xlabel('Time (days)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([-0.3 0.3])
title('u (m/s): linear fit, detrend, n2grid')


figure(9)
clf;set(gcf,'Color','w','Position', [41 277 1151 654]);
subplot(2,2,1)
pcolor(time_temp(plot_tidx),depth_n2,buoy1obs(plot_tidx,:)');shading flat;colorbar;
hold on;
contour(time_temp(plot_tidx),depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
hold off;
xlabel('Time (days)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
title('Buoyancy (m/s^2)')
colormap(WhiteBlueGreenYellowRed(0))
clim([-0.01 0.01]/4)

subplot(2,2,2)
pcolor(time_temp(plot_tidx),depth_reconst_n,n2_1obs(plot_tidx,:)');shading flat;colorbar;
hold on;
contour(time_temp(plot_tidx),depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
hold off;
xlabel('Time (days)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
title('N^2 (1/s^2)')
colormap(redblue)
clim([-1 1]/1e4/4)

subplot(2,2,3)
pcolor(time_temp(plot_tidx),depth_n2,buoy1fit(plot_tidx,:)');shading flat;colorbar;
hold on;
contour(time_temp(plot_tidx),depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
hold off;
xlabel('Time (days)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
title('Buoyancy (m/s^2): linear fit')
colormap(WhiteBlueGreenYellowRed(0))
clim([-0.01 0.01]/4)

subplot(2,2,4)
pcolor(time_temp(plot_tidx),depth_reconst_n,n2_1fit(plot_tidx,:)');shading flat;colorbar;
hold on;
contour(time_temp(plot_tidx),depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
hold off;
xlabel('Time (days)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
title('N^2 (1/s^2): linear fit')
colormap(redblue)
clim([-1 1]/1e4/4)


figure(1);
clf;set(gcf,'Color','w','Position',[114 662 1188*1.5 289]);
subplot(1,3,1)
pcolor(time_temp(plot_tidx),depth_temp,temp(plot_tidx,:)');shading flat;colorbar;
hold on;
contour(time_temp(plot_tidx),depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
hold off;
xlabel('Time (days)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([meanT-1 meanT+2])
title('Temperature (^oC)')

subplot(1,3,2)
pcolor(time_temp(plot_tidx),depth_temp,salt(plot_tidx,:)');shading flat;colorbar;
hold on;
contour(time_temp(plot_tidx),depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
hold off;
xlabel('Time (days)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([meanS-0.1 meanS+0.15])
title('Estimated salinity (psu)')

subplot(1,3,3)
pcolor(time_temp(plot_tidx),depth_n2,N2(plot_tidx,:)');shading flat;colorbar;
hold on;
contour(time_temp(plot_tidx),depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
hold off;
xlabel('Time (days)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
title('N^2 (1/s^2)')
clim([-1 1]/1e5)
colormap(cmocean('balance'));

