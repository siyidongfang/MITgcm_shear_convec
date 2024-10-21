clear all;close all
fontsize = 20;

fg1=figure(1);
clf;set(gcf,'Color','w','Position',[98 329 1000 800])

tiledlay = tiledlayout(2,2);
nexttile
load('fig4/Ri_flat.mat')
plot(1./Ri_min,shear_calc_Ri*1e3,'LineWidth',2)
hold on;
load('fig4/Ri_topo4_new.mat')
plot(1./Ri_min,shear_calc_Ri*1e3,'LineWidth',2)
grid on;grid minor;
xlim([0 10])

set(gca,'FontSize',fontsize)
legend('Flat bottom ($\theta=0^\circ$)','Sloping bottom ($\theta=4^\circ$)',...
    'Position', [0.5677 0.7325 0.2952 0.1],'interpreter','latex','Fontsize',fontsize+4);
legend('boxoff');
ylabel('Background shear $\Lambda$ (10$^{-3}\,$s $^{-1}$)','interpreter','latex');
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
title('${R_i}_\mathrm{min}^{-1}$ \textit{vs.} Background Shear','interpreter','latex','FontSize',fontsize+3);



nexttile(3)
load('../instability_km/exps_varyingN/N1e-3output','grow_rw','shear_all')
load('fig4/Ri_flat.mat')
[max_grow ~]=max(grow_rw,[],2);
for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km(i) = Ri_min(b(i));
end
plot(1./Ri_km,max_grow,'LineWidth',2);
hold on;
clear Ri_km max_grow shear_all shear_calc_Ri
load('../instability_km/exps_new/topo4_nu0_large_shearoutput.mat','grow_rw','shear_all')
load('fig4/Ri_topo4_new.mat')
[max_grow ~]=max(grow_rw,[],2);
for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km(i) = Ri_min(b(i));
end
plot(1./Ri_km,max_grow,'LineWidth',2);
grid on;grid minor;
xlim([0 10])
ylim([0 0.4])
set(gca,'FontSize',fontsize)
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
ylabel('(hour$^{-1}$)','interpreter','latex');
title('${R_i}_\mathrm{min}^{-1}$ \textit{vs.} Growth Rate','interpreter','latex','FontSize',fontsize+3);


nexttile(4)
load('../instability_km/exps_varyingN/N1e-3output','grow_rw','shear_all')
[max_grow ~]=max(grow_rw,[],2);
plot(shear_all*1e3,max_grow,'LineWidth',2);
hold on;
load('../instability_km/exps_new/topo4_nu0_large_shearoutput.mat','grow_rw','shear_all')
[max_grow ~]=max(grow_rw,[],2);
plot(shear_all*1e3,max_grow,'LineWidth',2);
grid on;grid minor;
xlim([0 3.33])
ylim([0 0.4])
set(gca,'FontSize',fontsize)
xlabel('Background shear $\Lambda$ (10$^{-3}\,$s $^{-1}$)','interpreter','latex');
ylabel('(hour$^{-1}$)','interpreter','latex');
title('Background Shear \textit{vs.} Growth Rate','interpreter','latex','FontSize',fontsize+3);



tiledlay.TileSpacing = 'compact';
tiledlay.Padding = 'compact';
AddLetters2Plots(fg1,'FontSize',fontsize+5,'FontWeight','normal')

print('-dpng','-r300',['fig_supp_new/figS_RivsShear.png']);

