
clear;
expdir = 'exps_Radko'

shear_Ri1 = 0.0031625;
shear_Ri5 = 0.001415;

shear = shear_Ri1;
m0_all = [0:0.015:6];
kx_all = [-0.1:0.0005:0.1];

% shear = shear_Ri5
% m0_all = [0:0.01:4];
% kx_all = [-0.5:0.0025:0.5];

growth=zeros(length(kx_all),length(m0_all));
for m=1:length(m0_all)
    m0=m0_all(m);
    load([expdir '/growth_shear' num2str(shear*1e3,3) '_m0' num2str(m0) '.mat'],'grow')
    growth(:,m) = grow; 
end


figure(1)
set(gcf,'Color','w')
pcolor(kx_all*0.01,m0_all*0.01,growth'/3600*1000);shading flat;colorbar;
colormap(jet);
% clim([0 0.015])
clim([0 6e-3])
xlim([-5 5]*1e-4)

set(gca,'Fontsize',20)
xlabel('Dimensionless $k = k^\star d$  (d=0.01 m)','Interpreter','latex')
ylabel('Dimensionless $m = m^\star d$  (d=0.01 m)','Interpreter','latex')
% title({'Dimensionless growth rate $\lambda = \lambda^\star \tau,\ (\tau=10^3\, s)$','$\overline{R_i}=1,\,\mathrm{dimensionless\ } \omega=\omega^\star\tau=0.1$','$N^2 = 1\times10^{-5}\,s^{-2},\ \Lambda = 0.0032\,s^{-1}$, tidal period = 17.45 hours'},'Interpreter','latex')

title('$\overline{R_i}=1,\,\mathrm{dimensionless\ } \omega=0.1$' ,'Interpreter','latex')

%%
% clear;
% 
% load('/Users/ysi/MITgcm_shear_convec/instability_km/output_Ri1/growth_topo0_NOdiff.mat')
% figure(1)
% semilogx(rw_all,growth)
% % plot(atand(1./rw_all),growth)
% grid on;grid minor;
% 
% 
% [maxgrowth,rwidx] = max(growth)
% rw_maxgrowth = rw_all(rwidx);
% load('/Users/ysi/MITgcm_shear_convec/instability_km/output_Ri1/topo0/mz1kx0.014962shear0.0031623_NOdiff.mat')
% plot_timeseires


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