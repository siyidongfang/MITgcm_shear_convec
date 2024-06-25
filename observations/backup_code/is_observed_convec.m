
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
% nStart = 1031-1;
% nEnd = nStart+24*2*4+3;
nStart = 263-1;
nEnd = nStart+24*2*4+3;
tidx = nStart:nEnd;
zidx = 4:18;
% uselect = u_reconstruct(zidx,tidx)';
uselect = uu(zidx,tidx)';
% wselect = ww_tilde(zidx,tidx)';

depth_u = depth(zidx);

meanU = mean(uselect,'omitnan');

time_u = time(tidx)/3600; %%% in hours
time_u = time_u-time_u(1);
time_u = time_u/24; %%% in days

Uzavg = mean((uselect),2);
figure(24)
clf;
plot(time_u,Uzavg/100);
hold on;
plot(time_u,shear_linear(tidx));

figure(25)
plot(meanU,depth_u,'LineWidth',2);axis ij;
grid on;
grid minor
set(gca,'Fontsize',fontsize);
set(gcf,'Color','w')
ylabel('Depth (m)')
title('Mean flow (m/s)')

uselect = uselect-mean(uselect);




%%
% Temperature data of two tidal cycles:
% temp = ncread('mavs2_20210718_level1.nc','__xarray_dataarray_variable__');
% depth_temp = ncread('mavs2_20210718_level1.nc','depth');
% time_temp = double(ncread('mavs2_20210718_level1.nc','time'));

addpath moorings_gridded_thermistors/level1/mavs2/
temp = ncread('mavs2_20210710.nc','__xarray_dataarray_variable__');
depth_temp = ncread('mavs2_20210710.nc','depth');
time_temp = double(ncread('mavs2_20210710.nc','time'));


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

uselect_n2grid = interp2(DD_u,TT_u,uselect,DD_n2,TT_n2,'linear');
% wselect_n2grid = interp2(DD_u,TT_u,wselect,DD_n2,TT_n2,'linear');

%%% time-mean, depth-averaged N^2 + re-constructed velocity
% adv1 = wselect_n2grid*meanN2+uselect_n2grid*meanN2*sind(topo);
adv1 = uselect_n2grid*meanN2*sind(topo);

%%% time-mean, depth-varying N^2 + re-constructed velocity
% adv2 = wselect_n2grid.*N2_zavg+uselect_n2grid.*N2_zavg*sind(topo);
% adv2 = uselect_n2grid.*N2_zavg*sind(topo);

%%% time-dependent, depth-varying N^2 + re-constructed velocity
% adv3 = wselect_n2grid.*N2+uselect_n2grid.*N2*sind(topo);
% adv3 = uselect_n2grid.*N2*sind(topo);

dt = 86400*(time_temp(2)-time_temp(1));

buoy1 =  - cumsum(adv1*dt,1);
% buoy2 =  - cumsum(adv2*dt,1);
% buoy3 =  - cumsum(adv3*dt,1);

% n20 = 0.5*(N2(1,1:end-1)+N2(1,2:end));

% n20 = mean(n20);
% n20 = meanN2;

Nt = 2*(24*3600+50*60);
avgN = 60;
% n20idx = [1:avgN Nt/4:Nt/4+avgN Nt/2:Nt/2+avgN Nt/4*3:Nt/4*3+avgN ...
%     Nt/8:Nt/8+avgN Nt/8*3:Nt/8*3+avgN Nt/8*5:Nt/8*5+avgN Nt/8*7*3:Nt/8*7+avgN];
% n20idx = [1:avgN Nt/4:Nt/4+avgN Nt/2:Nt/2+avgN Nt/4*3:Nt/4*3+avgN];
% n20 = mean(N2(n20idx,:),'all')

% n20z = mean(N2(n20idx,:),1);

[n20tz, depth_reconst_n] = gsw_Nsquared(0.5*(SA(1,1:end-1)+SA(1,2:end)),...
    0.5*(CT(1,1:end-1)+CT(1,2:end)),depth_n2,lat_MAVS2);

% buoy1 = smoothdata(buoy1,2,"gaussian");
% buoy2 = smoothdata(buoy2,2,"gaussian");
% buoy3 = smoothdata(buoy3,2,"gaussian");

n2_1 = meanN2*cosd(topo) + diff(buoy1,1,2)./diff(-depth_n2);
% n2_2 = n20 + diff(buoy2,1,2)./diff(-depth_n2);
% n2_3 = n20tz + diff(buoy3,1,2)./diff(-depth_n2);

plot_observed_convec;
