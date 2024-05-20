
load('MITgcm_growth_linearShear.mat')
figure(1)
clf;
set(gcf,'Color','w','Position',[55 396 956 404])
plot(shear_MITgcm,growth_MITgcm,'LineWidth',2)
hold on;
% plot(shear_MITgcm(4:13),growth_MITgcm2(4:13),'--','LineWidth',2)
set(gca,'Fontsize',20)
grid on;grid minor;
xlabel('Shear (1/s)')
ylabel('(1/hour)')
title('Growth rate (1/hour)')

load('../instability/GrowthRate_exps_linear.mat')
plot(shear_Floquet,growth_Floquet,':','LineWidth',3)

load('../instability_km/growth_experiments_flat_nu2e-4.mat')
plot(shear_km,growth_km,'-.','LineWidth',3)

legend('MITgcm, linear shear(\nu=\kappa=2\times10^{-4} m^2/s)',...
    'Floquet: periodic in x (\nu=\kappa=2\times10^{-4}m^2/s)',...
    'Floquet: double periodic (\nu=\kappa=2\times10^{-4} m^2/s)','Position', [0.1437 0.6645 0.4174 0.2438])

xlim([0 2.2]*1e-3)
