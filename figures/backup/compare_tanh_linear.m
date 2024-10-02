
clear;close all;
addpath ../analysis/colormaps/
fontsize = 17;
load_colors;

figure(1)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 950 900]);

%--- flat bottom
load('../instability_km/exps_new/topo0_nu0_output.mat')
load('fig4/Ri_flat.mat')

load('fig4/MIT_tanh_topo0_kv5e-6_tt.mat')
grow_gcm_tanh = growth_MITgcm;
shear_gcm_tanh = shear_MITgcm;
for i=1:length(shear_gcm_tanh)
    [a(i) b(i)] = min(abs(shear_MITgcm(i)-shear_calc_Ri));
    Ri_gcm_tanh(i) = Ri_min(b(i));
end
clear growth_MITgcm shear_MITgcm



load('fig4/MIT_flat_kv5e-6_tt.mat');
for i=1:length(shear_MITgcm)
    [a(i) b(i)] = min(abs(shear_MITgcm(i)-shear_calc_Ri));
    Ri_gcm(i) = Ri_min(b(i));
end

for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km(i) = Ri_min(b(i));
end

load('fig4/flat_km_kv2e-4.mat','max_grow','shear_all')
for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km_diff(i) = Ri_min(b(i));
end

load('../instability_km/exps_varyingN/N1e-3output','grow_rw','shear_all')
max_grow_km = max(grow_rw,[],2);
for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km(i) = Ri_min(b(i));
end

ax2 = subplot('position',[.57 .74 0.4 0.225]);
annotation('textbox',[0.538+0.02 0.996 0.15 0.01],'String','b','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
scatter(1./Ri_gcm,growth_MITgcm,36,blue,'LineWidth',2);
grid on;grid minor;
hold on;
scatter(1./Ri_gcm_tanh,grow_gcm_tanh,50,black,'LineWidth',2);
plot(1./Ri_km_diff,max_grow,'-','LineWidth',2,'Color',green);
plot(1./Ri_km,max_grow_km,'--','LineWidth',2,'Color',gray);
ylabel('(hour$^{-1}$)','interpreter','latex');
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
l1 = legend('MITgcm: linear shear', 'MITgcm: tanh shear','Linear theory: viscous', 'Linear theory: inviscid','Position',[0.5734 0.8749 0.2033 0.0877],'interpreter','latex');
set(gca,'Fontsize',fontsize);
xlim([0 8])
% ylim([-1e-3 0.4])
title('Growth rate (flat bottom)','interpreter','latex','Fontsize',fontsize+5);
box on;



