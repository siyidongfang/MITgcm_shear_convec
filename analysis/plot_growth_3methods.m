
% load('MITgcm_growth_linearShear.mat')
load('MITgcm_growth_tanhShear.mat')
figure(1)
clf;set(gcf,'Color','w','Position',[55 396 956 404])

Shear_all = (0:0.1:2)*1e-3;
N2 = 1e-6;
Ri_all = 1./(N2./(Shear_all).^2);
% xaxisvalue = Ri_all;
xaxisvalue = Shear_all;

plot(xaxisvalue,growth_MITgcm,'LineWidth',2)
hold on;
xlabel('Shear (1/s)')
 % xlabel('Ri^{-1}')

set(gca,'Fontsize',20)
grid on;grid minor;
ylabel('(1/hour)')
title('Growth rate (1/hour)')

load('../instability/GrowthRate_exps_tanh_BottomCenter_dz2.mat')
% plot(xaxisvalue,growth_Floquet,':','LineWidth',3)

growth_crop = max(GrowthRate_Floquet(4:end,:));
plot(shear_Floquet,growth_crop,'--','LineWidth',3)

% load('../instability_km/growth_experiments_flat_nu2e-4.mat')
% lam_x = 2*pi./kx_all;lam_z = 2*pi./m0_all;
% xstart = 7;zstart = 8;
% for s=1:length(shear_km)
%     growth_crop_km(s) = max(grow_smk(s,zstart:end,xstart:end),[],"all");
% end
% plot(xaxisvalue,growth_km,'-.','LineWidth',3)
% plot(xaxisvalue,growth_crop_km,'*','LineWidth',3)

legend('MITgcm, linear shear(\nu=\kappa=2\times10^{-4} m^2/s)',...
    'Floquet: periodic in x (\nu=\kappa=2\times10^{-4}m^2/s)',...
    'Floquet: double periodic (\nu=\kappa=2\times10^{-4} m^2/s)','Position', [0.1437 0.6645 0.4174 0.2438])

% xlim([0 2.2]*1e-3)
