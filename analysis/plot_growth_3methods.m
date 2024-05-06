
load('MITgcm_growth_hires_topo4.mat')
figure(1)
clf;
set(gcf,'Color','w','Position',[55 416 1080 384])
plot(shear_MITgcm,growth_MITgcm,'LineWidth',2)
hold on;
plot(shear_MITgcm,growth_MITgcm2,'--','LineWidth',2)
set(gca,'Fontsize',20)
grid on;grid minor;
xlabel('Shear (1/s)')
ylabel('(1/hour)')
title('Growth rate (1/hour)')

load('../instability/GrowthRate_Floquet_topo4.mat')
plot(shear_Floquet,growth_Floquet,':','LineWidth',3)
plot(shear_Floquet,growth_Floquet,':','LineWidth',3)

legend('MITgcm','MITgcm','Floquet: finite shear layer','Floquet: infinite shear layer','Position', [0.6056 0.2337 0.2389 0.2435])

