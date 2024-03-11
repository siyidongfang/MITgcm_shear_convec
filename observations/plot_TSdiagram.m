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
% load CTD.mat
% 
% lat_1to7 = mean(lat_all(1:7));
% lon_1to7 = mean(lon_all(1:7));
% 
% pt_1to7 = pt_all(1:823,1:7);
% CT_1to7 = CT_all(1:823,1:7);
% psal_1to7 = psal_all(1:823,1:7);
% SA_1to7 = SA_all(1:823,1:7);
% pp_1to7 = P_all(1:823,1:7);
% 
% psal_1to7 = psal_1to7(:)';
% pt_1to7 = pt_1to7(:)';
% pp_1to7 = pp_1to7(:)';
% %%% make T-S diagram
% p_ref = 0;
% pot_dens_contours = 20;
% make_TS_plot_gsw (psal_1to7,pt_1to7,pp_1to7,p_ref,pot_dens_contours,lat_1to7,lon_1to7)

load CTD_stations.mat

for nf = 1:15
    % nf = 1

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
    pot_dens_contours = 15;
    depths = pressure;
    [p1(nf) p2(nf)] = make_TS_plot_gsw (psal',pot_temp',depths',p_ref,pot_dens_contours,lat,lon);


    % print('-dpng','-r200',['CTD/fig_' num2str(nf) '.png']);

end

    %%% Linear fit T-S relation at MAVS2

    meanp1 = mean(p1(1:11));
    meanp2 = mean(p2(1:11));

    meanp1*5+meanp2;

