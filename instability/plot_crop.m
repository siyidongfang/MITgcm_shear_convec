fontsize = 20;
% load('GrowthRate_exps_tanh_ZeroCenter_dz2.mat')
load('GrowthRate_exps_tanh_BottomCenter_dz1.mat')
% load('GrowthRate_exps_linear_dz0.5.mat')

figure(3)
set(gcf,'color','w')
plot(shear_Floquet,growth_Floquet,'LineWidth',2);
grid on;grid minor;set(gca,'Fontsize',fontsize);
xlabel('Shear (1/s)')
ylabel('Growth rate (1/hour)')

figure(4)
set(gcf,'color','w')
pcolor(shear_Floquet,log10(lambda_Floquet),GrowthRate_Floquet);shading flat;colorbar;
clim([-0.3 0.3]);colormap(redblue)
grid on;grid minor;set(gca,'Fontsize',fontsize);
xlabel('Shear (1/s)')
title('Growth rate (1/hour)')
ylabel('log_{10}(\lambda_x) (m)')

crop_grow = GrowthRate_Floquet(4:end,:);
max_grow = max(crop_grow);
figure(5)
set(gcf,'color','w')
plot(shear_Floquet,max_grow,'LineWidth',2);
grid on;grid minor;set(gca,'Fontsize',fontsize);
xlabel('Shear (1/s)')
ylabel('Growth rate (1/hour)')