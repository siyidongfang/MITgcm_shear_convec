%%%
%%% calc_mean_dbdz.m
%%%
%%% Calculate the mean dbdz from observations
clear; close all;

addpath /Users/ysi/Software/gsw_matlab_v3_06_11/
addpath /Users/ysi/Software/gsw_matlab_v3_06_11/library/
addpath /Users/ysi/Software/gsw_matlab_v3_06_11/thermodynamics_from_t/


% Load the thermistor data

% %%% MAVS 1
% fpath = 'moorings_gridded_thermistors/level1/mavs1/mavs1_';
% Ndate = [801:2:829 902];
% lon_MAVS = -11.861534; 
% lat_MAVS = 54.198556; 
% load('MAVS1_velocity.mat','depth');
% u_zidx = 5:18;%%% For MAVS1
% depth=depth(u_zidx);
% n2_zidx = 9:44;%%% For MAVS1


%%% MAVS 2
fpath = 'moorings_gridded_thermistors/level1/mavs2/mavs2_';
Ndate = [706:2:720];
lon_MAVS = -11.843511;
lat_MAVS = 54.183718;
load('MAVS2_velocity.mat','depth');
u_zidx = 4:18;%%% For MAVS2
depth=depth(u_zidx);
n2_zidx =6:50;%%% For MAVS2


% idx = [n2_zidx(1) n2_zidx(end)];
idx = n2_zidx;

temp = [];
time_temp = [];
for n=1:length(Ndate)
    %%% MAVS 1
    % temp = [temp;ncread([fpath '20210' num2str(Ndate(n)) '.nc'],'t')];

    %%% MAVS 2
    temp = [temp;ncread([fpath '20210' num2str(Ndate(n)) '.nc'],'__xarray_dataarray_variable__')];
    
    time_temp = [time_temp;double(ncread([fpath '20210' num2str(Ndate(n)) '.nc'],'time'))];
end
depth_temp = ncread([fpath '20210' num2str(Ndate(1)) '.nc'],'depth');

temp=temp(:,idx);
depth_temp = depth_temp(idx);

%%% Calculate local db/dz
o1 = 0.078611022276845;
o2 = 34.638914089176580;
salt = o1*temp+o2;

[SA, ~] = gsw_SA_from_SP(salt,depth_temp,lon_MAVS,lat_MAVS);
CT = gsw_CT_from_pt(SA,temp);

N2 = NaN.*zeros(length(SA),length(depth_temp)-1);

smooth_SA = NaN*SA;
smooth_CT = NaN*CT;
smooth_N2 = NaN*N2;

Nsmooth = 900;

parfor kk=1:length(idx)
    smooth_SA(:,kk) = smoothdata(SA(:,kk),'gaussian',Nsmooth);
    smooth_CT(:,kk) = smoothdata(CT(:,kk),'gaussian',Nsmooth);
end

% rho = gsw_rho(SA,CT,depth_temp);
% smooth_rho = gsw_rho(smooth_SA,smooth_CT,depth_temp);
% g = 9.81;
% rho0 = 1000;
% N2_zavg = -g/rho0*(rho(:,2)-rho(:,1))./(-(depth_temp(2)-depth_temp(1)));
% smooth_N2_zavg = -g/rho0*(smooth_rho(:,2)-smooth_rho(:,1))./(-(depth_temp(2)-depth_temp(1)));

parfor ii = 1:length(SA)
    [N2(ii,:), ~] = gsw_Nsquared(SA(ii,:)',CT(ii,:)',depth_temp,lat_MAVS);
    [smooth_N2(ii,:), ~] = gsw_Nsquared(smooth_SA(ii,:)',smooth_CT(ii,:)',depth_temp,lat_MAVS);
end

N2_zavg = mean(N2,2,'omitnan');
smooth_N2_zavg = mean(smooth_N2,2,'omitnan');

save('MAVS2_N2.mat','N2_zavg','time_temp','smooth_N2_zavg')









