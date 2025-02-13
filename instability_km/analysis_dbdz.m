
% clear;
% load('output_Ri1_Nsq1e-6_topo4/growth_topo4_NOdiff.mat')
% figure(1)
% semilogx(rw_all,growth)
% 
% [maxgrow rwidx] = max(growth)
% 
% rw_maxgrow = rw_all(rwidx);
% kx = mz*rw_maxgrow;
% 
% load('/Users/ysi/MITgcm_shear_convec/instability_km/output_Ri1_Nsq1e-6_topo4/mz0.01kx0.0015849shear0.00097.mat')
% plot_timeseires

clear;

shear_all = [0:0.1:1.8]*1e-3;
NS = length(shear_all);

expdir = 'output/topo4_Nsq1e-6';

% for ns = 1:NS
for ns =19
    shear = shear_all(ns);
    % load(['output/topo4_Nsq1e-6/growth_shear' num2str(shear*1e3,3) '.mat'])
    % [max_growth(ns) rw_idx(ns)] = max(grow);
    % rw_mg(ns) = rw_all(rw_idx(ns));
    load([expdir '/growth_shear' num2str(shear*1e3,3) '_analysis.mat'])
    plot_timeseires
    print('-djpeg','-r150',[expdir '/growth_shear' num2str(shear*1e3,3) '_analysis.jpeg']);
    % grow
    kx/m0
end



