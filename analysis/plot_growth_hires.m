
addpath functions/

clear;
% close all;
list_exps;
for ne = 1:21
    expname = EXPNAME{ne};
    loadexp;
    load([expdir expname '/RMSE.mat'])
    growth(ne) = pp(1);
    load([expdir expname '/RMSE2.mat'])
    growth2(ne) = pp(1);
    if(growth(ne)>0.1)
        growth2(ne) = growth(ne);
    end
end

shear_all = [0:0.1:2]*1e-3;
figure(1)
clf;
set(gcf,'Color','w')
plot(shear_all,growth,'LineWidth',2)
hold on;
plot(shear_all,growth2,'--','LineWidth',2)
set(gca,'Fontsize',20)
grid on;grid minor;
xlabel('Shear (1/s)')
ylabel('(1/hour)')
title('Growth rate (1/hour)')

growth_MITgcm = growth;
growth_MITgcm2 = growth2;
shear_MITgcm = shear_all;

save('MITgcm_growth_hires_topo4.mat','growth_MITgcm','growth_MITgcm2','shear_MITgcm')
