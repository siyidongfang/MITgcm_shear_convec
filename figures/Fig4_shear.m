
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
% load('../instability_km/exps_new/topo0_nu0_output.mat')
load('../instability_km/exps_varyingN/N1e-3output.mat')
load('fig4/Ri_flat.mat')

[max_grow_nodiff r_idx]=max(grow_rw,[],2);
for i=1:length(shear_all)
    r_mostunstable(i) = 1./rw_all(r_idx(i));
end
% r_mostunstable(r_mostunstable>17)=NaN;

for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km(i) = Ri_min(b(i));
end

%%
ax1 = subplot('position',[.065 .74 0.375 0.225]);
annotation('textbox',[0.028+0.02 0.996 0.15 0.01],'String','A','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(shear_all,1./rw_all,grow_rw');
hold on;
scatter(shear_all,r_mostunstable,7,brown,'filled','LineWidth',1);
shading interp;
set(gca,'Fontsize',fontsize);
colormap(WhiteBlueGreenYellowRed(0))
h1 = colorbar;
set(h1,'Position',[0.45 0.7456+0.01   0.008    0.2]);
set(get(h1,'Title'),'String',{'$\ \ \ \ (\mathrm{hour}^{-1})$'},'interpreter','latex','FontSize',fontsize);
clim([0 0.4]);
xlabel('Tidal shear $\Lambda$ (s$^{-1}$)','interpreter','latex');
ylabel('Perturbation aspect ratio $m_0/k_0$','interpreter','latex');
title('Growth rate (flat bottom)','interpreter','latex','Fontsize',fontsize+5);
ylim([0 30])
xlim([0 2.1]*1e-3)
set(gca,'TickDir','out');
grid on;grid minor;

%%%%%%%%%%Start loading
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


load('fig4/flat_km_kv2e-4.mat','max_grow','shear_all')
shear_diff = shear_all;
for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km_diff(i) = Ri_min(b(i));
end

load('../instability_km/exps_varyingN/N1e-3output','grow_rw','shear_all')
shear_km = shear_all;
max_grow_km = max(grow_rw,[],2);
for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km(i) = Ri_min(b(i));
end

ax2 = subplot('position',[.57 .74 0.4 0.225]);
annotation('textbox',[0.538+0.02 0.996 0.15 0.01],'String','B','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
scatter(shear_MITgcm,growth_MITgcm,36,blue,'LineWidth',2);
grid on;grid minor;
hold on;
scatter(shear_gcm_tanh,grow_gcm_tanh,50,black,'LineWidth',2);
plot(shear_diff,max_grow,'-','LineWidth',2,'Color',green);
plot(shear_km,max_grow_km,'--','LineWidth',2,'Color',darkgray);
ylabel('(hour$^{-1}$)','interpreter','latex');
xlabel('Tidal shear $\Lambda$ (s$^{-1}$)','interpreter','latex');
l1 = legend('MITgcm: linear shear', 'MITgcm: tanh shear','Linear theory: viscous', 'Linear theory: inviscid','Position',[0.5734 0.8749 0.2033 0.0877],'interpreter','latex');
set(gca,'Fontsize',fontsize);
% xlim([0 8])
% ylim([-1e-3 0.4])
title('Growth rate (flat bottom)','interpreter','latex','Fontsize',fontsize+5);
box on;
xlim([0 2.1]*1e-3)

%%%%%%%%%%End loading

% ax2 = subplot('position',[.57 .74 0.4 0.225]);
% annotation('textbox',[0.538+0.02 0.996 0.15 0.01],'String','b','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
% scatter(1./Ri_gcm,growth_MITgcm,36,blue,'LineWidth',2);
% grid on;grid minor;
% hold on;
% plot(1./Ri_km_diff,max_grow,'-','LineWidth',2,'Color',green);
% % plot(1./Ri_km,max_grow_nodiff,'--','LineWidth',2,'Color',yellow);
% ylabel('(hour$^{-1}$)','interpreter','latex');
% xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
% l1 = legend('MITgcm','Theory','Position',[0.58 0.91 0.1871 0.0458],'interpreter','latex');
% set(gca,'Fontsize',fontsize);
% xlim([0 5.2])
% ylim([-1e-3 0.4])
% title('Growth rate (flat bottom)','interpreter','latex','Fontsize',fontsize+5);
% box on;





%%
clear growth_MITgcm max_grow_rw Ri_gcm Ri_km shear_MITgcm shear_all shear_calc_Ri
%--- topo = 4 degrees
% load('../instability_km/exps_new/topo4_nu0_output.mat')
load('../instability_km/exps_new/topo4_nu0_large_shearoutput.mat')
shear_km = shear_all;
load('fig4/Ri_topo4.mat')

lam_x_real = lam_x_real/6;
max_grow_rw = max(grow_rw,[],2);

load('fig4/MIT_zeroCenter_tanh_topo4_kv5e-6_tt.mat')
grow_gcm_tanh = growth_MITgcm;
shear_gcm_tanh = shear_MITgcm;
for i=1:length(shear_gcm_tanh)
    [a(i) b(i)] = min(abs(shear_MITgcm(i)-shear_calc_Ri));
    Ri_gcm_tanh(i) = Ri_min(b(i));
end
clear growth_MITgcm shear_MITgcm


load('fig4/MIT_topo4_kv5e-6_tt.mat')
for i=1:length(shear_MITgcm)
    if(shear_MITgcm(i)>shear_convec)
        shear_MITgcm(i)=NaN;
        growth_MITgcm(i)=NaN;
        Ri_gcm(i)=NaN;
    else
        [a(i) b(i)] = min(abs(shear_MITgcm(i)-shear_calc_Ri));
        Ri_gcm(i) = Ri_min(b(i));
    end
end

for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km(i) = Ri_min(b(i));
end

clear r_mostunstable
[max_grow r_idx]=max(grow_rw,[],2);
for i=1:length(shear_all)
    r_mostunstable(i) = 1./rw_all(r_idx(i));
end
r_mostunstable(r_mostunstable>12)=NaN;



clear max_grow shear_all Ri_km_diff 
load('fig4/topo4_km_kv2e-4.mat','max_grow','shear_all')
shear_diff = shear_all;
for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km_diff(i) = Ri_min(b(i));
end





ax3 = subplot('position',[.06 .40 0.375 0.225]);
annotation('textbox',[0.028+0.02 0.65 0.15 0.01],'String','C','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
rw_idx = 1:Nrw;
pcolor(shear_km,1./rw_all,grow_rw');
hold on;
scatter(shear_km,r_mostunstable,7,brown,'filled','LineWidth',1);
plot(shear_km,zeros(1,length(Ri_km)),'k--','LineWidth',0.5);
set(gca,'Fontsize',fontsize);
shading interp;
h2 = colorbar;
set(h2,'Position',[0.45 0.4056+0.01   0.008    0.2]);
set(get(h2,'Title'),'String',{'$\ \ \ \ (\mathrm{hour}^{-1})$'},'interpreter','latex','FontSize',fontsize);
clim([0 0.4]);
xlabel('Tidal shear $\Lambda$ (s$^{-1}$)','interpreter','latex');
ylabel('Perturbation aspect ratio $m_0/k_0$','interpreter','latex');
title('Growth rate (sloping bottom)','interpreter','latex','Fontsize',fontsize+5);
ylim([-21 14])
% xlim([0 8])
xlim([0 2.1]*1e-3)
set(gca,'TickDir','out');
grid on;grid minor;

max_grow_km = max(grow_rw,[],2);

ax4 = subplot('position',[.57 .40 0.4 0.225]);
annotation('textbox',[0.538+0.02 0.65 0.15 0.01],'String','D','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
scatter(shear_MITgcm,growth_MITgcm,36,blue,'LineWidth',2);
grid on;grid minor;
hold on;
% scatter(1./Ri_gcm_tanh,grow_gcm_tanh,50,black,'LineWidth',2);
plot(shear_diff,max_grow,'-','LineWidth',2,'Color',green);
plot(shear_km,max_grow_km,'--','LineWidth',2,'Color',darkgray);
ylabel('(hour$^{-1}$)','interpreter','latex');
xlabel('Tidal shear $\Lambda$ (s$^{-1}$)','interpreter','latex');
% l4 = legend('MITgcm: linear shear','MITgcm: tanh shear','Linear theory: viscous','Linear theory: inviscid','Position',[0.58 0.57 0.1871 0.0458],'interpreter','latex');
l4 = legend('MITgcm: linear shear','Linear theory: viscous','Linear theory: inviscid','Position',[0.5751 0.5517 0.2033 0.0668],'interpreter','latex');
set(gca,'Fontsize',fontsize);
% xlim([0 5.2])
% xlim([0 8])
% ylim([-1e-3 0.6])
xlim([0 2.1]*1e-3)
ylim([0 0.4])
title('Growth rate (sloping bottom)','interpreter','latex','Fontsize',fontsize+5);
box on;



print('-dpng','-r300','fig4/fig4_shear.png');

