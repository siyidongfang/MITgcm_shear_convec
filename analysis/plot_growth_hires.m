
addpath functions/

clear;close all;
list_exps;

for ne = [1:10 21]
    expname = EXPNAME{ne};
    loadexp;
    load([expdir expname '/RMSE_tt.mat'])
    growth(ne) = pp(1);
end

growth=squeeze(growth([1:10 21]));

shear_MITgcm = [0:0.1:0.9 2.0]*1e-3;
% shear_MITgcm = [0:0.1:2.0]*1e-3;
 % shear_MITgcm = [0:0.1:1.8]*1e-3;
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

save('../figures/fig4/MIT_flat.mat','growth_MITgcm','shear_MITgcm')


%%

% load('MITgcm_growth_hires_flat.mat')
figure(1)
clf
set(gcf,'Color','w')
plot(shear_MITgcm,growth_MITgcm,'LineWidth',2)
hold on;
set(gca,'Fontsize',20)
grid on;grid minor;
xlabel('Shear (1/s)')
ylabel('(1/hour)')
title('Growth rate (1/hour)')

% load('MITgcm_growth_hires_topo4.mat')
% plot(shear_MITgcm,growth_MITgcm,'--','LineWidth',2)
% 
% legend('Flat bottom', 'topography = 4 degrees')




%%

