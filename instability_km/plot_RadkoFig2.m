
clear;

topo = 0;
load(['output_Nsquare1e-5/growth_topo' num2str(topo) '.mat'])

figure(1)
pcolor(growth);shading flat;colorbar;colormap(jet);clim([0 0.3])
% topo = 0;
% mz_all = [0:0.25:6];
% mz_plot = repmat(mz_all,[41 1]);

% load(['output_Nsquare1e-5/growth_topo' num2str(topo) '_mz_' num2str(mz) '.mat'])

% kx_plot = rw_all'.*mz_all;
% 
% figure(1)
% pcolor(log10(kx_plot));colorbar;shading flat;
% 
% idx = 11;
% shear_all(idx)
% growth_Ri1 = growth(idx,:);
% semilogx(kx,growth_Ri1);