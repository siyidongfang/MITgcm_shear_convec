
clear;close all;
addpath ../instability/
addpath ../analysis/colormaps/

load_colors;
showplot = true;
fontsize = 20;

%%% Velocity data:
load('MAVS2_LinearShear.mat');

time_adcp = double(ncread('MAVS2_24606.nc','time')) + 33.000323*1e6; % microseconds since 2021-07-07 06:00:33.000323
Time_adcp = datetime(2021,7,7) + hours(round(time_adcp/1e6/3600)+6);
Time_adcp.Format = 'yyyy-MM-dd';

% Find time indices in the ADCP data of two tidal cycles:
nStart = 1031;
nEnd = 1222;
% nStart = 263;
% nEnd = 454;
tidx = nStart:nEnd;
zidx = 4:18;
uselect = u_reconstruct(zidx,tidx)';
% uselect = uu_tilde(zidx,tidx)';

time_u = time(tidx)/3600; %%% in hours
time_u = time_u-time_u(1);
time_u = time_u/24; %%% in days

depth_u = depth(zidx);

%%% Temperature data of two tidal cycles:
temp = ncread('mavs2_20210718_level1.nc','__xarray_dataarray_variable__');
depth_temp = ncread('mavs2_20210718_level1.nc','depth');
time_temp = double(ncread('mavs2_20210718_level1.nc','time'));
time_temp = time_temp'/3600;%%% in hours
time_temp = time_temp/24;%%% in days

% addpath moorings_gridded_thermistors/level1/mavs2/
% temp = ncread('mavs2_20210710.nc','__xarray_dataarray_variable__');
% depth_temp = ncread('mavs2_20210710.nc','depth');
% time_temp = double(ncread('mavs2_20210710.nc','time'));

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

N2_zavg = mean(N2,'omitnan');
meanN2 = mean(N2_zavg);


%%% Interpolate low-resolution velocity data onto the grid of high-res
%%% N2 grid
[DD_u,TT_u] = meshgrid(depth_u,time_u);
[DD_n2,TT_n2] = meshgrid(depth_n2,time_temp);

uselect_n2grid = interp2(DD_u,TT_u,uselect,DD_n2,TT_n2,'linear');

%%% time-mean, depth-averaged N^2 + re-constructed velocity
adv1 = uselect_n2grid*meanN2*sind(topo);

%%% time-mean, depth-varying N^2 + re-constructed velocity
adv2 = uselect_n2grid.*N2_zavg*sind(topo);

%%% time-dependent, depth-varying N^2 + re-constructed velocity
adv3 = uselect_n2grid.*N2*sind(topo);

dt = 86400*(time_temp(2)-time_temp(1));

buoy1 =  - cumsum(adv1*dt,1);
buoy2 =  - cumsum(adv2*dt,1);
buoy3 =  - cumsum(adv3*dt,1);

% n20 = 0.5*(N2(1,1:end-1)+N2(1,2:end));
[n20, depth_reconst_n] = gsw_Nsquared(0.5*(SA(1,1:end-1)+SA(1,2:end)),...
    0.5*(CT(1,1:end-1)+CT(1,2:end)),depth_n2,lat_MAVS2);

n2_1 = n20 + diff(buoy1,1,2)./diff(-depth_n2);
n2_2 = n20 + diff(buoy2,1,2)./diff(-depth_n2);
n2_3 = n20 + diff(buoy3,1,2)./diff(-depth_n2);

% plot_observed_convec;
