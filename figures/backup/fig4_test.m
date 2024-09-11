

% clear;close all;
addpath ../analysis/colormaps/
fontsize = 17;
load_colors;

clear growth_MITgcm max_grow_rw Ri_gcm Ri_km shear_MITgcm shear_all shear_calc_Ri
%--- topo = 4 degrees
load('backup/MIT_topo4_hires_kv1e-5.mat')
% load('../instability_km/exps_new/topo4_nu0_output.mat')
load('fig4/Ri_topo4.mat')

 crop_limit = 32000;

rw_idx_crop=find(lam_x_real<=crop_limit);
max_grow_rw = max(grow_rw(:,rw_idx_crop),[],2);

for i=1:length(shear_MITgcm)
    [a(i) b(i)] = min(abs(shear_MITgcm(i)-shear_calc_Ri));
    Ri_gcm(i) = Ri_min(b(i));
end

for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km(i) = Ri_min(b(i));
end

figure(1)
% plot(1./Ri_km,max_grow_rw,'LineWidth',2,'Color',black);
grid on;grid minor;
hold on;
plot(1./Ri_gcm,growth_MITgcm,'LineWidth',2,'Color',red);
ylabel('(hour$^{-1}$)','interpreter','latex');
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
set(gca,'Fontsize',fontsize);
xlim([0 4])
ylim([-1e-3 0.35])
title('Growth rate (sloping bottom)','interpreter','latex','Fontsize',fontsize+5);


l4 = legend('Theory','MITgcm, 2e-4','MITgcm, 1e-4','MITgcm, 5e-5','MITgcm, 1e-5','interpreter','latex');







