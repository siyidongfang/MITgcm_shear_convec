
addpath functions/

clear;
% close all;
list_exps;
for ne = 1:19
    expname = EXPNAME{ne};
    loadexp;
    load([expdir expname '/RMSE.mat'])
    growth(ne) = pp(1);
end

shear_MITgcm = [0:0.1:1.8]*1e-3;
figure(1)
% clf;
set(gcf,'Color','w')
plot(shear_MITgcm,growth,'LineWidth',2)
hold on;
set(gca,'Fontsize',20)
grid on;grid minor;
xlabel('Shear (1/s)')
ylabel('(1/hour)')
title('Growth rate (1/hour)')

growth_MITgcm = growth;

save('MITgcm_growth_hires_topo4.mat','growth_MITgcm','shear_MITgcm')
