
% clear;close all;
addpath ../analysis/colormaps/
fontsize = 17;
load_colors;

figure(1)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 950 900]);

%--- flat bottom
% load('fig4/MIT_exps_flat_test.mat')
% load('fig4/MITgcm_growth_hires_flat.mat')
load('fig4/MIT_flat_kv8e-5_tt.mat')

load('../instability_km/exps_new/topo0_nu0_output.mat')
load('fig4/Ri_flat.mat')
load('../instability/products/GrowthRate_exps_linear_dz0.5.mat')

lam_x_real=lam_x_real/6;
crop_limit = 40000/6;
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

load('fig4/Floquet_km_flat_8e-5.mat')
% load('fig4/flat_km_kv2e-4.mat','max_grow','shear_all')
for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km_diff(i) = Ri_min(b(i));
end

ax1 = subplot('position',[.065 .74 0.375 0.225]);
annotation('textbox',[0.028+0.02 0.996 0.15 0.01],'String','a','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
rw_idx = 1:Nrw;
pcolor(1./Ri_km,(lam_x_real(rw_idx))/1000,grow_rw(:,rw_idx)');
shading interp;
hold on;
% plot(1./Ri_km,crop_limit/1000*ones(1,length(Ri_km)),'k--','LineWidth',1);
hold off;
set(gca,'Fontsize',fontsize);
% colormap(flip(bone))
% colormap(cmocean('tempo'))
% colormap(cmocean('rain'))
colormap(WhiteBlueGreenYellowRed(0))
h1 = colorbar;
set(h1,'Position',[0.45 0.7456+0.01   0.008    0.2]);
set(get(h1,'Title'),'String',{'$\ \ \ \ (\mathrm{hour}^{-1})$'},'interpreter','latex','FontSize',fontsize);
clim([0 0.35]);
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
ylabel('Horizontal wavelength (km)','interpreter','latex');
title('Growth rate (flat bottom)','interpreter','latex','Fontsize',fontsize+5);
ylim([0 33]/6)
xlim([0 4])

load('fig4/Floquet_km_flat.mat')
ax2 = subplot('position',[.57 .74 0.4 0.225]);
annotation('textbox',[0.538+0.02 0.996 0.15 0.01],'String','b','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
plot(1./Ri_gcm,growth_MITgcm,'LineWidth',2,'Color',blue);
grid on;grid minor;
hold on;
plot(1./Ri_gcm,max_grow_floquet,'LineWidth',2,'Color',black);
plot(1./Ri_gcm,max_grow_floquet_exclude,'LineWidth',2,'Color',red);
% plot(1./Ri_km,max_grow_rw,'LineWidth',2,'Color',black);
plot(1./Ri_km_diff,max_grow,'--','LineWidth',2,'Color',green);
% % load('../instability_km/exps_new/topo0_nu2e-4_lores_output.mat')
% % crop_limit = 3000;
% % rw_idx=find(lam_x_real<=crop_limit);
% % [max_grow,I]= max(grow_rw(:,rw_idx),[],2);
% % plot(1./Ri_gcm,max_grow,'LineWidth',2,'Color',red);
% grow_250m = GrowthRate_Floquet(20,:);
% grow_320m = GrowthRate_Floquet(18,:);
% grow_400m = GrowthRate_Floquet(16,:);
% plot(1./Ri_gcm,grow_250m,'-.','LineWidth',1.5)
% plot(1./Ri_gcm,grow_320m,'-.','LineWidth',1.5)
% plot(1./Ri_gcm,grow_400m,'-.','LineWidth',1.5)
% plot(1./Ri_gcm,growth_Floquet,':','LineWidth',2)
ylabel('(hour$^{-1}$)','interpreter','latex');
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
l1 = legend('MITgcm','Theory: Eigenvalues','Theory: cutoff Lx=380m','Theory: Time Advancement','Position',[0.5807 0.8896 0.2427 0.0668],'interpreter','latex');
% l1 = legend('MITgcm','Theory','Position',[0.58 0.9140 0.1010 0.0445],'interpreter','latex');
% l1 = legend('MITgcm','Theory','Discretized, Lx = 250m','Discretized, Lx = 320m','Discretized, Lx = 400m','Discretized, Lx =300$\sim$700m','interpreter','latex');
set(gca,'Fontsize',fontsize);
xlim([0 4])
ylim([-1e-3 0.35])
title('Growth rate (flat bottom)','interpreter','latex','Fontsize',fontsize+5);



%%
clear growth_MITgcm max_grow_rw Ri_gcm Ri_km shear_MITgcm shear_all shear_calc_Ri
%--- topo = 4 degrees
% load('fig4/MIT_topo4_hires_kv1e-4.mat')
% load('fig4/MITgcm_growth_hires_topo4.mat')
load('fig4/MIT_exps_topo4.mat')
load('../instability_km/exps_new/topo4_nu0_output.mat')
load('fig4/Ri_topo4.mat')

lam_x_real = lam_x_real/6;
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

clear max_grow shear_all Ri_km_diff
load('fig4/topo4_km_kv2e-4.mat','max_grow','shear_all')
for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km_diff(i) = Ri_min(b(i));
end

ax3 = subplot('position',[.06 .40 0.375 0.225]);
annotation('textbox',[0.028+0.02 0.65 0.15 0.01],'String','c','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
rw_idx = 1:Nrw;
pcolor(1./Ri_km,(lam_x_real(rw_idx))/1000,grow_rw(:,rw_idx)');
hold on;
% plot(1./Ri_km,crop_limit/1000*ones(1,length(Ri_km)),'k--','LineWidth',1);
hold off;
set(gca,'Fontsize',fontsize);
shading interp;
h2 = colorbar;
set(h2,'Position',[0.45 0.4056+0.01   0.008    0.2]);
set(get(h2,'Title'),'String',{'$\ \ \ \ (\mathrm{hour}^{-1})$'},'interpreter','latex','FontSize',fontsize);
clim([0 0.35]);
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
ylabel('Across-slope wavelength (km)','interpreter','latex');
title('Growth rate (sloping bottom)','interpreter','latex','Fontsize',fontsize+5);
ylim([0 33]/6)
xlim([0 4])


load('fig4/Floquet_km_topo4.mat')
ax4 = subplot('position',[.57 .40 0.4 0.225]);
annotation('textbox',[0.538+0.02 0.65 0.15 0.01],'String','d','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
plot(1./Ri_gcm,growth_MITgcm,'LineWidth',2,'Color',blue);
grid on;grid minor;
hold on;
plot(1./Ri_gcm,max_grow_floquet,'LineWidth',2,'Color',black);
plot(1./Ri_gcm,max_grow_floquet_exclude,'LineWidth',2,'Color',red);
plot(1./Ri_km_diff,max_grow,'--','LineWidth',2,'Color',green);
% plot(1./Ri_km,max_grow_rw,'--','LineWidth',2,'Color',brown);
ylabel('(hour$^{-1}$)','interpreter','latex');
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
l4 = legend('MITgcm','Theory: Eigenvalues','Theory: cutoff Lx=380m','Theory: Time Advancement','Position',[0.5807 0.5501 0.2427 0.0668],'interpreter','latex');
% l4 = legend('MITgcm','Theory: inviscid','Theory: $\kappa=\nu=2e-4$','Position',[0.58 0.5712 0.1010 0.0445],'interpreter','latex');
% l4 = legend('MITgcm','Theory','Position',[0.58 0.5712 0.1010 0.0445],'interpreter','latex');
set(gca,'Fontsize',fontsize);
xlim([0 4])
ylim([-1e-3 0.35])
title('Growth rate (sloping bottom)','interpreter','latex','Fontsize',fontsize+5);




load('fig4/fig4_km.mat')
tidx = nt_percycle*5+1:nt_percycle*10;
%- Calculate the buoyancy budget
uB0x = -re_uuu(tidx)*N^2*ss;
wB0z = -re_www(tidx)*N^2*cs;
wBz  =  re_www(tidx)*shear/omega*N^2*ss.*st(tidx);
% diffusion = 0*
dbdt = [0 (re_buoy(3:end)-re_buoy(1:end-2))/dt/2 0];
dbdt = dbdt(tidx);

%- Normalization
uB0x = uB0x/max(abs(dbdt));
wB0z = wB0z/max(abs(dbdt));
wBz = wBz/max(abs(dbdt));
dbdt = dbdt/max(abs(dbdt));
tt = tt(tidx);


%--- plot vertical velocity and buoyancy perturbation
ax5 = subplot('position',[.06 .06 0.4 0.225]);
annotation('textbox',[0.028+0.02 0.31 0.15 0.01],'String','e','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
plot(tt/43200,re_buoy(tidx)./max(abs(re_buoy(tidx))),'LineWidth',2,'Color',blue);
hold on;
plot(tt/43200,re_www(tidx)./max(abs(re_www(tidx))),'LineWidth',2,'Color',black);
grid on;grid minor;
l5 = legend('Buoyancy perturbation','Vertical velocity','interpreter','latex','Position',[0.07 0.2265 0.2063 0.0458]);
xlabel('Time (tidal cycles)','interpreter','latex')
set(gca,'Fontsize',fontsize);
title('Normalized perturbations','interpreter','latex','Fontsize',fontsize+5);


ax6 = subplot('position',[.57 .06 0.4 0.225]);
annotation('textbox',[0.538+0.02 0.31 0.15 0.01],'String','f','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
lres = plot(tt/43200,dbdt-uB0x-wB0z-wBz,'-','LineWidth',3,'Color',boxcolor);
hold on;
luB0x = plot(tt/43200,uB0x,'LineWidth',2,'Color',brown);
lwBz = plot(tt/43200,wBz,'LineWidth',2,'Color',yellow);
lwB0z = plot(tt/43200,wB0z,'LineWidth',2,'Color',blue);
ldbdt = plot(tt/43200,dbdt,'--','LineWidth',1,'Color',black);
grid on;grid minor;
xlabel('Time (tidal cycles)','interpreter','latex')
set(gca,'Fontsize',fontsize);
title('Normalized buoyancy budget','interpreter','latex','Fontsize',fontsize+5);
ylim([-1 1])
l61 = legend([lwB0z,lwBz,luB0x], ...
    '$w^\prime\partial_z B_0$','$w^\prime\partial_z B$','$u^\prime\partial_x B_0$',...
    'interpreter','latex','Position',[0.58 0.2104 0.0889 0.0668]);
ah=axes('position',get(ax6,'position'),'visible','off');
set(gca,'Fontsize',fontsize);
l62 = legend(ah,[ldbdt lres], ...
    '$\partial_\tau b^\prime = \partial_t b^\prime + U\partial_x b^\prime$','Residual',...
    'interpreter','latex','Position', [0.58 0.0825 0.1626 0.0458]');
legend('boxoff') 


% print('-dpng','-r300','fig4/fig4_matlab.png');