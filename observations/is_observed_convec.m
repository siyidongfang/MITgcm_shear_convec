
clear;close all;
addpath ../instability/
addpath ../analysis/colormaps/

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
% nStart = 263;
% nEnd = 454;
tidx = nStart:nEnd;
zidx = 4:18;
uselect = uu(zidx,tidx)';
vselect = vv(zidx,tidx)';
wselect = ww(zidx,tidx)';
time_uw = time(tidx)/3600; %%% in hours
time_uw = time_uw-time_uw(1);

depth_uw = depth(zidx);

if(size(find(isnan(uselect)))~=0)
    uselect(3,8)=0.5*(uselect(2,8)+uselect(4,8));
    wselect(3,8)=0.5*(wselect(2,8)+wselect(4,8));
end
% 
% windowsize = round(400./mean(diff(depth_uw)));
% uselect = smoothdata(uselect,2,'gaussian',windowsize,'omitnan');

%%% Temperature data of two tidal cycles:
temp = ncread('mavs2_20210718_level1.nc','__xarray_dataarray_variable__');
depth_temp = ncread('mavs2_20210718_level1.nc','depth');
time_temp = double(ncread('mavs2_20210718_level1.nc','time'));
% addpath moorings_gridded_thermistors/level1/mavs2/
% temp = ncread('mavs2_20210710.nc','__xarray_dataarray_variable__');
% depth_temp = ncread('mavs2_20210710.nc','depth');
% time_temp = double(ncread('mavs2_20210710.nc','time'));

% temp = temp(21600:end,:);
% time_temp = time_temp(21600:end);

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

%%% Add the contribution of salinity to buoyancy:
load('CTD/CTD.mat','SA_all','P_all','SA_mean15','p_mid15')
% salinity = SA_all(:,9);
% depth_salinity = P_all(:,9);
sBeta = 1e-3;
salinity = SA_mean15;
depth_salinity = P_all(:,9);

%%% Interpolate salinity to buoyancy grid
salt_uw_stagger1 = interp1(depth_salinity(590:711),salinity(590:711),depth_uw_stagger');
salt_uw_stagger1 = salt_uw_stagger1';
figure(20)
plot(salinity,depth_salinity);axis ij
hold on;
plot(salt_uw_stagger1,depth_uw_stagger)

salt_uw_stagger = repmat(salt_uw_stagger1,[length(temp_uw_stagger) 1]);

N2_uwgrid = gravity*tAlpha.*diff(temp_uw_stagger,1,2)./repmat(diff(-depth_uw_stagger),[length(temp_uw_stagger) 1]);

N2_uwgrid = N2_uwgrid - gravity*sBeta.*diff(salt_uw_stagger,1,2)./repmat(diff(-depth_uw_stagger),[length(salt_uw_stagger) 1]);

Ttavg_uw_stagger= mean(temp_uw_stagger);
N2avg_uwgrid = mean(N2_uwgrid);


% ttselect = 40:192;
ttselect = 1:length(uselect);


%%% time-mean, depth-averaged N^2 + prescribed velocities
uprescribe = zeros(size(uselect));
adv0 = uprescribe*meanN2*sind(topo);
 
%%% time-mean, depth-averaged N^2 + observed velocities
adv1 = uselect(ttselect,:)*meanN2*sind(topo)+wselect(ttselect,:)*meanN2*cosd(topo);

%%% time-mean, depth-varying N^2 + observed velocities
adv2 = uselect(ttselect,:).*N2avg_uwgrid*sind(topo)+wselect(ttselect,:).*N2avg_uwgrid*cosd(topo);

%%% observed N^2 + observed velocities
adv3 = uselect(ttselect,:).*N2_uwgrid(ttselect,:)*sind(topo)+wselect(ttselect,:).*N2_uwgrid(ttselect,:)*cosd(topo);

% %%% time-mean, depth-averaged N^2 + observed velocities
% adv1 = wselect(ttselect,:)*meanN2*cosd(topo);
% 
% %%% time-mean, depth-varying N^2 + observed velocities
% adv2 = wselect(ttselect,:).*N2avg_uwgrid*cosd(topo);
% 
% %%% observed N^2 + observed velocities
% adv3 = wselect(ttselect,:).*N2_uwgrid(ttselect,:)*cosd(topo);


% %%% time-mean, depth-averaged N^2 + observed velocities
% adv1 = uselect(ttselect,:)*meanN2*sind(topo);
% 
% %%% time-mean, depth-varying N^2 + observed velocities
% adv2 = uselect(ttselect,:).*N2avg_uwgrid*sind(topo);
% 
% %%% observed N^2 + observed velocities
% adv3 = uselect(ttselect,:).*N2_uwgrid(ttselect,:)*sind(topo);

% b0 = 0;
b01 = gravity*tAlpha*0.5*(temp_uw_stagger(1,1:end-1)+temp_uw_stagger(1,2:end));
b0 = gravity*(tAlpha*0.5*(temp_uw_stagger(1,1:end-1)+temp_uw_stagger(1,2:end))...
    -sBeta*0.5*(salt_uw_stagger(1,1:end-1)+salt_uw_stagger(1,2:end)));

dt = 3600*(time_uw(2)-time_uw(1));

buoy0 = b0 - cumsum(adv0*dt,1);
buoy1 = b0 - cumsum(adv1*dt,1);
buoy2 = b0 - cumsum(adv2*dt,1);
buoy3 = b0 - cumsum(adv3*dt,1);

n2_0 = diff(buoy0,1,2)./diff(-depth_uw);
n2_1 = diff(buoy1,1,2)./diff(-depth_uw);
n2_2 = diff(buoy2,1,2)./diff(-depth_uw);
n2_3 = diff(buoy3,1,2)./diff(-depth_uw);
depth_n2grid_recons = 0.5*(depth_uw(1:end-1)+depth_uw(2:end));


plot_observed_convec;
