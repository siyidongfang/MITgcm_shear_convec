%%%
%%% plot_TSdiagram.m
%%%
%%% plot the T-S diagram using CTD data

%%% load the data
clear;close all;
addpath /Users/ysi/Software/gsw_matlab_v3_06_11/thermodynamics_from_t/;
addpath /Users/ysi/Software/gsw_matlab_v3_06_11/library/;
addpath /Users/ysi/Software/gsw_matlab_v3_06_11/;
addpath ../analysis/colormaps/
addpath CTD
load CTD_stations.mat

for nf = 1:15

Nz = length(CTD.Stn(nf).pressure);
pressure = CTD.Stn(nf).pressure;
insitu_temp = CTD.Stn(nf).insitu_temp;
psal = CTD.Stn(nf).psal;
lat = CTD.Stn(nf).lat;
lon = CTD.Stn(nf).lon;

%%% Calculate absolute salinity
[SA, in_ocean] = gsw_SA_from_SP(psal,pressure,lon,lat);

%%% Calculate potential temperature
p_ref = 0;
pot_temp = gsw_pt_from_t(SA,insitu_temp,pressure,p_ref);


%%% make T-S diagram
pot_dens_contours = 20;
depths = pressure;
make_TS_plot_gsw (psal',pot_temp',depths',p_ref,pot_dens_contours,lat,lon)

end
% load_colors
