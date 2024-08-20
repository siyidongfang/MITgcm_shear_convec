
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
% load('../instability/products/GrowthRate_exps_linear_dz0.5.mat')
% load('fig4/MIT_flat_kv8e-5_tt.mat');
% load('fig4/MIT_flat_noCori_tt.mat');
load('fig4/MIT_flat_kv5e-6_tt.mat');
% load('fig4/Floquet_km_flat_8e-5.mat')

% lam_x_real=lam_x_real/6;
% max_grow_rw = max(grow_rw,[],2);

for i=1:length(shear_MITgcm)
    [a(i) b(i)] = min(abs(shear_MITgcm(i)-shear_calc_Ri));
    Ri_gcm(i) = Ri_min(b(i));
end

% for i=1:length(shear_floquet)
%     [a(i) b(i)] = min(abs(shear_floquet(i)-shear_calc_Ri));
%     Ri_floquet(i) = Ri_min(b(i));
% end

for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km(i) = Ri_min(b(i));
end



% load('fig4/flat_km_kv2e-4.mat','max_grow','shear_all')
% for i=1:length(shear_all)
%     [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
%     Ri_km_diff(i) = Ri_min(b(i));
% end

load('fig5/fig5_topo0_noDiff_eig.mat')
for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_eig(i) = Ri_min(b(i));
end
max_grow_eig = max(grow_all);
rw_eig = lam_z_all./lam_x_all;

ax1 = subplot('position',[.065 .74 0.375 0.225]);
annotation('textbox',[0.028+0.02 0.996 0.15 0.01],'String','a','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
% pcolor(1./Ri_km,lam_x_real/1000,grow_rw');
% pcolor(1./Ri_km,1./rw_all,grow_rw');
pcolor(1./Ri_eig,1./rw_eig,grow_all);
hold on;

% figure(2)
% inverseRi = (omega/N*1./rw_all).^2;
% plot(inverseRi,1./rw_all);
xxxplot = 1./Ri_eig;
xxxplot(xxxplot<2)=NaN;
yyyplot = N/omega.*sqrt(xxxplot-2);
plot(xxxplot,yyyplot,'LineWidth',2,'Color',black);
xxxplot = 1./Ri_eig;
xxxplot(xxxplot<2.5)=NaN;
yyyplot = N/omega.*sqrt(xxxplot-2.5);
plot(xxxplot,yyyplot,'LineWidth',2,'Color',black);
xxxplot = 1./Ri_eig;
xxxplot(xxxplot<3)=NaN;
yyyplot = N/omega.*sqrt(xxxplot-3);
plot(xxxplot,yyyplot,'LineWidth',2,'Color',black);
xxxplot = 1./Ri_eig;
xxxplot(xxxplot<3.5)=NaN;
yyyplot = N/omega.*sqrt(xxxplot-3.5);
plot(xxxplot,yyyplot,'LineWidth',2,'Color',black);
xxxplot = 1./Ri_eig;
yyyplot = N/omega.*sqrt(xxxplot+0.8);
plot(xxxplot,yyyplot,'LineWidth',2,'Color',black);
% xxxplot = 1./Ri_km;
% yyyplot = N/omega.*sqrt(xxxplot)+1;
% plot(xxxplot,yyyplot,'LineWidth',2,'Color',black);
% xlim([0 4])
% ylim([0 20])

% figure(1)

% set(gca, 'XScale', 'log')
shading interp;
set(gca,'Fontsize',fontsize);
colormap(WhiteBlueGreenYellowRed(0))
h1 = colorbar;
set(h1,'Position',[0.45 0.7456+0.01   0.008    0.2]);
set(get(h1,'Title'),'String',{'$\ \ \ \ (\mathrm{hour}^{-1})$'},'interpreter','latex','FontSize',fontsize);
clim([0 0.4]);
% xlabel('Minimum Richardson number ${R_i}_\mathrm{min}$','interpreter','latex');
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
% ylabel('Horizontal wavelength (km)','interpreter','latex');
ylabel('Perturbation aspect ratio $m_0/k_0$','interpreter','latex');
title('Growth rate (flat bottom)','interpreter','latex','Fontsize',fontsize+5);
% ylim([0 33]/6)
ylim([-18 18])

xlim([0  5.2])
% xticks([0.1 0.2 0.3 0.5 0.75 1 2 3 4])




%%
ax2 = subplot('position',[.57 .74 0.4 0.225]);
annotation('textbox',[0.538+0.02 0.996 0.15 0.01],'String','b','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
scatter(1./Ri_gcm,growth_MITgcm,36,blue,'LineWidth',2);
grid on;grid minor;
hold on;
% scatter(Ri_floquet,max_grow_floquet,36,black,'LineWidth',2);
% plot(Ri_gcm,growth_MITgcm,'-','LineWidth',2,'Color',blue);
plot(1./Ri_eig,max_grow_eig,'-','LineWidth',2,'Color',black);
% plot(1./Ri_floquet,max_grow_floquet,'-','LineWidth',2,'Color',black);
% set(gca, 'XScale', 'log')
% scatter(shear_floquet,max_grow_floquet,'LineWidth',2,'Color',black);
% scatter(shear_MITgcm,growth_MITgcm,'LineWidth',2,'Color',blue);
% plot(1./Ri_km_diff,max_grow,'--','LineWidth',2,'Color',green);
ylabel('(hour$^{-1}$)','interpreter','latex');
% xlabel('Minimum Richardson number ${R_i}_\mathrm{min}$','interpreter','latex');
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
% l1 = legend('MITgcm','Theory: Eigenvalues','Theory: Time Advancement','Position',[0.5807 0.8896 0.2427 0.0668],'interpreter','latex');
l1 = legend('MITgcm','Theory: Eigenvalues','Position',[0.58 0.91 0.1871 0.0458],'interpreter','latex');
set(gca,'Fontsize',fontsize);
xlim([0 5.2])
ylim([-1e-3 0.4])
title('Growth rate (flat bottom)','interpreter','latex','Fontsize',fontsize+5);
box on;

%%
clear growth_MITgcm max_grow_rw Ri_gcm Ri_km shear_MITgcm shear_all shear_calc_Ri
%--- topo = 4 degrees
load('../instability_km/exps_new/topo4_nu0_output.mat')
load('fig4/Ri_topo4.mat')

lam_x_real = lam_x_real/6;
max_grow_rw = max(grow_rw,[],2);

% load('fig4/MIT_topo4_noCori_tt.mat')
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



clear max_grow shear_all Ri_km_diff 
load('fig4/topo4_km_kv2e-4.mat','max_grow','shear_all')
for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km_diff(i) = Ri_min(b(i));
end


load('fig5/fig5_topo4_noDiff_eig.mat')
for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_eig(i) = Ri_min(b(i));
end
max_grow_eig = max(grow_all);
rw_eig = lam_z_all./lam_x_all;

ax3 = subplot('position',[.06 .40 0.375 0.225]);
annotation('textbox',[0.028+0.02 0.65 0.15 0.01],'String','c','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
rw_idx = 1:Nrw;
% pcolor(1./Ri_km,(lam_x_real(rw_idx))/1000,grow_rw(:,rw_idx)');
% pcolor(1./Ri_km,1./rw_all,grow_rw');
pcolor(1./Ri_eig,1./rw_eig,grow_all);
% set(gca, 'XScale', 'log')
set(gca,'Fontsize',fontsize);
shading interp;
h2 = colorbar;
set(h2,'Position',[0.45 0.4056+0.01   0.008    0.2]);
set(get(h2,'Title'),'String',{'$\ \ \ \ (\mathrm{hour}^{-1})$'},'interpreter','latex','FontSize',fontsize);
clim([0 0.4]);
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
% ylabel('Across-slope wavelength (km)','interpreter','latex');
ylabel('Perturbation aspect ratio $m_0/k_0$','interpreter','latex');
title('Growth rate (sloping bottom)','interpreter','latex','Fontsize',fontsize+5);
ylim([-18 18])
xlim([0 5.2])


% clear Ri_floquet shear_floquet 
% load('fig4/Floquet_km_topo4_8e-5.mat')
% shear_floquet=shear_all;
% for i=1:length(shear_floquet)
%     [a(i) b(i)] = min(abs(shear_floquet(i)-shear_calc_Ri));
%     Ri_floquet(i) = Ri_min(b(i));
% end



ax4 = subplot('position',[.57 .40 0.4 0.225]);
annotation('textbox',[0.538+0.02 0.65 0.15 0.01],'String','d','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
scatter(1./Ri_gcm,growth_MITgcm,36,blue,'LineWidth',2);
grid on;grid minor;
hold on;
% scatter(1./Ri_floquet,max_grow_floquet,36,black,'LineWidth',2);
plot(1./Ri_eig,max_grow_eig,'-','LineWidth',2,'Color',black);
% plot(1./Ri_floquet,max_grow_floquet,'-','LineWidth',2,'Color',black);
% plot(1./Ri_gcm,growth_MITgcm,'-','LineWidth',2,'Color',blue);
% set(gca, 'XScale', 'log')
% scatter(shear_floquet,max_grow_floquet,'LineWidth',2,'Color',black);
% scatter(shear_MITgcm,growth_MITgcm,'LineWidth',2,'Color',blue);
% plot(1./Ri_km_diff,max_grow,'--','LineWidth',2,'Color',green);
ylabel('(hour$^{-1}$)','interpreter','latex');
% xlabel('Minimum Richardson number ${R_i}_\mathrm{min}$','interpreter','latex');
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
% l4 = legend('MITgcm','Theory: Eigenvalues','Theory: Time Advancement','Position',[0.5807 0.5501 0.2427 0.0668],'interpreter','latex');
l4 = legend('MITgcm','Theory: Eigenvalues','Position',[0.58 0.57 0.1871 0.0458],'interpreter','latex');
set(gca,'Fontsize',fontsize);
xlim([0 5.2])
ylim([-1e-3 0.4])
title('Growth rate (sloping bottom)','interpreter','latex','Fontsize',fontsize+5);
box on;


load('fig4/fig4_km.mat')
tidx = nt_percycle*5+1:nt_percycle*10;
%- Calculate the buoyancy budget
uB0x = -re_uuu(tidx)*N^2*ss;
wB0z = -re_www(tidx)*N^2*cs;
wBz  =  re_www(tidx)*shear/omega*N^2*ss.*st(tidx);
% diffusion = 0
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
ldbdt = plot(tt/43200,dbdt,'LineWidth',1,'Color',black);
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

