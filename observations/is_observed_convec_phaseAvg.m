
clear;close all;
addpath ../instability/
fontsize = 20;
load('MAVS2_Phase_Average.mat')

hab = flip(hab)';
N2 = flip(N2)';
temp = flip(temp)';

meanT = mean(temp,'all','omitnan');
T_tavg = mean(temp,'omitnan');
N2avg = mean(N2,'omitnan');

zmax = 1500;
depth_temp = zmax-hab(1,:);
time_temp = phase./2/pi*12.42;
depth_n2 = depth_temp;

tAlpha = 2e-4;
gravity = 9.81;

load('MAVS2_velocity.mat','uu','vv','ww','depth','time','topo');

% nStart = 263+22;
nStart = 1031-5;
nEnd = nStart+48;
tidx = nStart:nEnd;
zidx = 4:18;
uselect = uu(zidx,tidx)';
vselect = vv(zidx,tidx)';
wselect = ww(zidx,tidx)';
time_uw = time(tidx)/3600; %%% in hours
time_uw = time_uw-time_uw(1);

depth_uw = depth(zidx);

%%% Interpolate velocity and stratification on the same 
%%% vertical grid and time grid
% depth_uw2 = depth(min(zidx)-1:max(zidx)+1);
% depth_uw_stagger = 0.5*(depth_uw2(1:end-1)+depth_uw2(2:end));

[DD_temp,TT_temp] = meshgrid(depth_temp,time_temp);
% [DD_uw_stagger,TT_uw_stagger] = meshgrid(depth_uw_stagger,time_uw);
[DD_uw,TT_uw] = meshgrid(depth_uw,time_uw);

% temp_uw_stagger = interp2(DD_temp,TT_temp,temp,DD_uw_stagger,TT_uw_stagger);
% N2_uwgrid = gravity*tAlpha.*diff(temp_uw_stagger,1,2)./repmat(diff(-depth_uw_stagger),[length(temp_uw_stagger) 1]);
temp_uwgrid = interp2(DD_temp,TT_temp,temp,DD_uw,TT_uw);
N2_uwgrid = interp2(DD_temp,TT_temp,N2,DD_uw,TT_uw);

depth_uw_stagger = depth_uw;
temp_uw_stagger = temp_uwgrid;

Ttavg_uw= mean(temp_uwgrid,'omitnan');
N2avg_uwgrid = mean(N2_uwgrid,'omitnan');

meanN2 = mean(N2avg);

meanT = mean(temp,'all','omitnan');

ttselect = 1:length(uselect);

windowsize = round(400./mean(diff(depth_uw)));
uselect = smoothdata(uselect,2,'gaussian',windowsize,'omitnan');

%%
%%% time-mean, depth-averaged N^2 + observed velocities
adv1 = uselect*meanN2*sind(topo)+wselect*meanN2*cosd(topo);

%%% time-mean, depth-varying N^2 + observed velocities
adv2 = uselect.*N2avg_uwgrid*sind(topo)+wselect.*N2avg_uwgrid*cosd(topo);

%%% observed N^2 + observed velocities
adv3 = uselect.*N2_uwgrid*sind(topo)+wselect.*N2_uwgrid*cosd(topo);


% %%% time-mean, depth-averaged N^2 + observed velocities
% adv1 = wselect*meanN2*cosd(topo);
% 
% %%% time-mean, depth-varying N^2 + observed velocities
% adv2 = wselect.*N2avg_uwgrid*cosd(topo);
% 
% %%% observed N^2 + observed velocities
% adv3 = wselect.*N2_uwgrid*cosd(topo);


% %%% time-mean, depth-averaged N^2 + observed velocities
% adv1 = uselect*meanN2*sind(topo);
% 
% %%% time-mean, depth-varying N^2 + observed velocities
% adv2 = uselect.*N2avg_uwgrid*sind(topo);
% 
% %%% observed N^2 + observed velocities
% adv3 = uselect.*N2_uwgrid*sind(topo);

% b0 =0;

b0 = gravity*tAlpha*temp_uwgrid(2,:);
dt = 3600*(time_uw(2)-time_uw(1));

buoy1 = b0 - cumsum(adv1*dt,1,'omitnan');
buoy2 = b0 - cumsum(adv2*dt,1,'omitnan');
buoy3 = b0 - cumsum(adv3*dt,1,'omitnan');

n2_1 = diff(buoy1,1,2)./diff(-depth_uw);
n2_2 = diff(buoy2,1,2)./diff(-depth_uw);
n2_3 = diff(buoy3,1,2)./diff(-depth_uw);
depth_n2grid_recons = 0.5*(depth_uw(1:end-1)+depth_uw(2:end));

n2_1(:,1)=NaN;
n2_2(:,1)=NaN;
n2_3(:,1)=NaN;

plot_observed_convec;


