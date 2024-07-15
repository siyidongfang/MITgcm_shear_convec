
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
load('fig3/MITgcm_growth_hires_flat.mat')
load('../instability_km/exps_new/topo0_nu0_output.mat')
load('fig3/Ri_flat.mat')

crop_limit = 4000;

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

ax1 = subplot('position',[.065 .74 0.375 0.225]);
annotation('textbox',[0.028 0.996 0.15 0.01],'String','a','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
rw_idx = 1:Nrw;
pcolor(1./Ri_km,(lam_x_real(rw_idx))/1000,grow_rw(:,rw_idx)');
shading interp;
hold on;
plot(1./Ri_km,crop_limit/1000*ones(1,length(Ri_km)),'k--','LineWidth',1);
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
ylim([0 33])
xlim([0 4])


ax2 = subplot('position',[.57 .74 0.4 0.225]);
annotation('textbox',[0.538 0.996 0.15 0.01],'String','b','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
plot(1./Ri_gcm,growth_MITgcm,'LineWidth',2,'Color',blue);
grid on;grid minor;
hold on;
plot(1./Ri_km,max_grow_rw,'LineWidth',2,'Color',black);
ylabel('(hour$^{-1}$)','interpreter','latex');
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
l1 = legend('MITgcm','Theory','Position',[0.58 0.9140 0.1010 0.0445],'interpreter','latex');
set(gca,'Fontsize',fontsize);
xlim([0 4])
ylim([-1e-3 0.35])
title('Growth rate (flat bottom)','interpreter','latex','Fontsize',fontsize+5);




clear growth_MITgcm max_grow_rw Ri_gcm Ri_km shear_MITgcm shear_all shear_calc_Ri
%--- topo = 4 degrees
load('fig3/MITgcm_growth_hires_topo4.mat')
load('../instability_km/exps_new/topo4_nu0_output.mat')
load('fig3/Ri_topo4.mat')

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


ax3 = subplot('position',[.06 .40 0.375 0.225]);
annotation('textbox',[0.028 0.65 0.15 0.01],'String','c','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
rw_idx = 1:Nrw;
pcolor(1./Ri_km,(lam_x_real(rw_idx))/1000,grow_rw(:,rw_idx)');
hold on;
plot(1./Ri_km,crop_limit/1000*ones(1,length(Ri_km)),'k--','LineWidth',1);
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
ylim([0 33])
xlim([0 4])


ax4 = subplot('position',[.57 .40 0.4 0.225]);
annotation('textbox',[0.538 0.65 0.15 0.01],'String','d','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
plot(1./Ri_gcm,growth_MITgcm,'LineWidth',2,'Color',blue);
grid on;grid minor;
hold on;
plot(1./Ri_km,max_grow_rw,'LineWidth',2,'Color',black);
ylabel('(hour$^{-1}$)','interpreter','latex');
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
l4 = legend('MITgcm','Theory','Position',[0.58 0.5712 0.1010 0.0445],'interpreter','latex');
set(gca,'Fontsize',fontsize);
xlim([0 4])
ylim([-1e-3 0.35])
title('Growth rate (sloping bottom)','interpreter','latex','Fontsize',fontsize+5);




load('fig3/fig3_km.mat')
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
annotation('textbox',[0.028 0.31 0.15 0.01],'String','e','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
plot(tt/43200,re_buoy(tidx)./max(abs(re_buoy(tidx))),'LineWidth',2,'Color',blue);
hold on;
plot(tt/43200,re_www(tidx)./max(abs(re_www(tidx))),'LineWidth',2,'Color',black);
grid on;grid minor;
l5 = legend('Buoyancy perturbation','Vertical velocity','interpreter','latex','Position',[0.07 0.2265 0.2063 0.0458]);
xlabel('Time (tidal cycles)','interpreter','latex')
set(gca,'Fontsize',fontsize);
title('Normalized perturbations','interpreter','latex','Fontsize',fontsize+5);


ax6 = subplot('position',[.57 .06 0.4 0.225]);
annotation('textbox',[0.538 0.31 0.15 0.01],'String','f','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
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


print('-dpng','-r300','fig3/fig3_matlab.png');