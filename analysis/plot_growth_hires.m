
addpath functions/

clear;close all;
list_exps;

for ne = 1:14
    expname = EXPNAME{ne};
    loadexp;
    load([expdir expname '/RMSE_tt.mat'])

    growth(ne) = pp(1);
end

shear_MITgcm = [0.1 0.4:0.2:2.8]*1e-03; % topo0, tanh
% shear_MITgcm = [0.1 0.4 0.6 0.8 1.0:0.1:1.3 1.335 1.4:0.1:2.0 2.08]*1e-03;%% topo4
% shear_MITgcm = [0.1 0.4 0.6 0.8 1.0:0.1:1.4 1.455 1.5:0.1:2.2]*1e-03; %% flat

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

save('../figures/fig4/MIT_tanh_topo0_kv5e-6_tt.mat','growth_MITgcm','shear_MITgcm')


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

