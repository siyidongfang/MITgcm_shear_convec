
clear;

topo = 0;
% load(['output_Ri1/growth_topo' num2str(topo) '_test3.mat'])
load('/Users/ysi/MITgcm_shear_convec/instability_km/output_Ri1_Nsq1e-6_topo4/growth_topo4.mat')


figure(1)
set(gcf,'Color','w')
pcolor(kx_all*0.01,mz_all*0.01,growth/3600*1000);shading flat;colorbar;
colormap(jet);
% clim([0.001 0.015])
set(gca,'Fontsize',20)
xlabel('Dimensionless $k = k^\star d$  (d=0.01 m)','Interpreter','latex')
ylabel('Dimensionless $m = m^\star d$  (d=0.01 m)','Interpreter','latex')
title({'Dimensionless growth rate $\lambda = \lambda^\star \tau,\ (\tau=10^3\, s)$','$\overline{R_i}=1,\,\mathrm{dimensionless\ } \omega=\omega^\star\tau=0.1$','$N^2 = 1\times10^{-5}\,s^{-2},\ \Lambda = 0.0032\,s^{-1}$, tidal period = 17.45 hours'},'Interpreter','latex')



%%
clear;

load('/Users/ysi/MITgcm_shear_convec/instability_km/output_Ri1/growth_topo0_NOdiff.mat')
figure(1)
semilogx(rw_all,growth)
% plot(atand(1./rw_all),growth)
grid on;grid minor;


[maxgrowth,rwidx] = max(growth)
rw_maxgrowth = rw_all(rwidx);
load('/Users/ysi/MITgcm_shear_convec/instability_km/output_Ri1/topo0/mz1kx0.014962shear0.0031623_NOdiff.mat')
plot_timeseires


%%
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