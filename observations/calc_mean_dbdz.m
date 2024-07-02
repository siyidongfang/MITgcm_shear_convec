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

%%% MAVS 2
fpath = 'moorings_gridded_thermistors/level1/mavs2/mavs2_';
Ndate = [706:2:720];
lon_MAVS = -11.843511;
lat_MAVS = 54.183718;

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


%%% Calculate local db/dz
o1 = 0.078611022276845;
o2 = 34.638914089176580;
salt = o1*temp+o2;

[SA, ~] = gsw_SA_from_SP(salt,depth_temp,lon_MAVS,lat_MAVS);
CT = gsw_CT_from_pt(SA,temp);

N2 = NaN.*zeros(length(SA),length(depth_temp)-1);

clear temp salt 

parfor ii = 1:length(SA)
    [N2(ii,:), depth_n2] = gsw_Nsquared(SA(ii,:)',CT(ii,:)',depth_temp,lat_MAVS);
end

%%% Calculate depth-averaged db/dz
N2_zavg = mean(N2,2,'omitnan');
N2_zavg_withNAN = mean(N2,2);

save('MAVS2_N2.mat','N2_zavg','time_temp','N2_zavg_withNAN')









