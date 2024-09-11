
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
load('fig4/MIT_flat_kv5e-6_tt.mat');

for i=1:length(shear_MITgcm)
    [a(i) b(i)] = min(abs(shear_MITgcm(i)-shear_calc_Ri));
    Ri_gcm(i) = Ri_min(b(i));
end

for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km(i) = Ri_min(b(i));
end

[max_grow_nodiff r_idx]=max(grow_rw,[],2);
for i=1:length(shear_all)
    r_mostunstable(i) = 1./rw_all(r_idx(i));
end
% r_mostunstable(r_mostunstable>17)=NaN;

%%
load('fig4/flat_km_kv2e-4.mat','max_grow','shear_all')
for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km_diff(i) = Ri_min(b(i));
end


ax1 = subplot('position',[.065 .74 0.375 0.225]);
annotation('textbox',[0.028+0.02 0.996 0.15 0.01],'String','a','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(1./Ri_km,1./rw_all,grow_rw');
hold on;
scatter(1./Ri_km,r_mostunstable,8,brown,'filled','LineWidth',1);
shading interp;
set(gca,'Fontsize',fontsize);
colormap(WhiteBlueGreenYellowRed(0))
h1 = colorbar;
set(h1,'Position',[0.45 0.7456+0.01   0.008    0.2]);
set(get(h1,'Title'),'String',{'$\ \ \ \ (\mathrm{hour}^{-1})$'},'interpreter','latex','FontSize',fontsize);
clim([0 0.4]);
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
ylabel('Perturbation aspect ratio $m_0/k_0$','interpreter','latex');
title('Growth rate (flat bottom)','interpreter','latex','Fontsize',fontsize+5);
ylim([0 22])
xlim([0 5.2])


ax2 = subplot('position',[.57 .74 0.4 0.225]);
annotation('textbox',[0.538+0.02 0.996 0.15 0.01],'String','b','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
scatter(1./Ri_gcm,growth_MITgcm,36,blue,'LineWidth',2);
grid on;grid minor;
hold on;
plot(1./Ri_km_diff,max_grow,'-','LineWidth',2,'Color',green);
% plot(1./Ri_km,max_grow_nodiff,'--','LineWidth',2,'Color',yellow);
ylabel('(hour$^{-1}$)','interpreter','latex');
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
l1 = legend('MITgcm','Theory','Position',[0.58 0.91 0.1871 0.0458],'interpreter','latex');
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
for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km_diff(i) = Ri_min(b(i));
end


ax3 = subplot('position',[.06 .40 0.375 0.225]);
annotation('textbox',[0.028+0.02 0.65 0.15 0.01],'String','c','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
rw_idx = 1:Nrw;
pcolor(1./Ri_km,1./rw_all,grow_rw');
hold on;
scatter(1./Ri_km,r_mostunstable,8,brown,'filled','LineWidth',1);
set(gca,'Fontsize',fontsize);
shading interp;
h2 = colorbar;
set(h2,'Position',[0.45 0.4056+0.01   0.008    0.2]);
set(get(h2,'Title'),'String',{'$\ \ \ \ (\mathrm{hour}^{-1})$'},'interpreter','latex','FontSize',fontsize);
clim([0 0.4]);
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
ylabel('Perturbation aspect ratio $m_0/k_0$','interpreter','latex');
title('Growth rate (sloping bottom)','interpreter','latex','Fontsize',fontsize+5);
ylim([0 22])
xlim([0 5.2])


ax4 = subplot('position',[.57 .40 0.4 0.225]);
annotation('textbox',[0.538+0.02 0.65 0.15 0.01],'String','d','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
scatter(1./Ri_gcm,growth_MITgcm,36,blue,'LineWidth',2);
grid on;grid minor;
hold on;
plot(1./Ri_km_diff,max_grow,'-','LineWidth',2,'Color',green);
ylabel('(hour$^{-1}$)','interpreter','latex');
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
l4 = legend('MITgcm','Theory','Position',[0.58 0.57 0.1871 0.0458],'interpreter','latex');
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

%--- Compute the vorticity budget
re_zeta = real(zeta);
dzetadt = [0 (re_zeta(3:end)-re_zeta(1:end-2))/dt/2 0];
dzetadt = dzetadt(tidx);

bxcs = real(1i*kx.*buoy(tidx)*cs);
bzss = real(-1i*mz_t(tidx).*buoy(tidx)*ss);

%--- Normalization
bxcs = bxcs/max(abs(dzetadt));
bzss = bzss/max(abs(dzetadt));
dzetadt = dzetadt/max(abs(dzetadt));


ax5 = subplot('position',[.06 .06 0.4 0.225]);
annotation('textbox',[0.028+0.02 0.31 0.15 0.01],'String','e','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
lzres = plot(tt/43200,dzetadt-bxcs-bzss,'-','LineWidth',3,'Color',boxcolor);
hold on;
lbxcs = plot(tt/43200,bxcs,'LineWidth',2,'Color',blue);
lbzss = plot(tt/43200,bzss,'LineWidth',2,'Color',yellow);
ldzetadt = plot(tt/43200,dzetadt,'LineWidth',1,'Color',green);
grid on;grid minor;
xlabel('Time (tidal cycles)','interpreter','latex')
set(gca,'Fontsize',fontsize);
title('Normalized vorticity budget','interpreter','latex','Fontsize',fontsize+5);
ylim([-1 1])
l51 = legend([lbxcs,lbzss], ...
    '$\partial_x b^\prime \cos\theta$','-$\partial_z b^\prime \sin\theta$',...
    'interpreter','latex','Position',[0.075 0.225 0.0889 0.0668/4*3]);
ah=axes('position',get(ax5,'position'),'visible','off');
set(gca,'Fontsize',fontsize);
l52 = legend(ah,[ldzetadt lzres], ...
    '$\partial_\tau \zeta = \partial_t \zeta+ U\partial_x \zeta$','Residual',...
    'interpreter','latex','Position', [0.065 0.0825 0.1626 0.0458]');
legend('boxoff') 



ax6 = subplot('position',[.57 .06 0.4 0.225]);
annotation('textbox',[0.538+0.02 0.31 0.15 0.01],'String','f','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
lres = plot(tt/43200,dbdt-uB0x-wB0z-wBz,'-','LineWidth',3,'Color',boxcolor);
hold on;
luB0x = plot(tt/43200,uB0x,'LineWidth',2,'Color',RED1);
lwBz = plot(tt/43200,wBz,'LineWidth',2,'Color',yellow);
lwB0z = plot(tt/43200,wB0z,'LineWidth',2,'Color',blue);
ldbdt = plot(tt/43200,dbdt,'LineWidth',1,'Color',green);
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



figure(2)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 950 900]);

%--- plot vertical velocity and buoyancy perturbation
ax5 = subplot('position',[.06 .06 0.4 0.225]);
annotation('textbox',[0.028+0.02 0.31 0.15 0.01],'String','e','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
plot(tt/43200,re_buoy(tidx)./max(abs(re_buoy(tidx))),'LineWidth',2,'Color',blue);
hold on;
plot(tt/43200,re_www(tidx)./max(abs(re_www(tidx))),'LineWidth',2,'Color',black);
plot(tt/43200,re_zeta(tidx)./max(abs(re_zeta(tidx))),'LineWidth',2,'Color',green);
grid on;grid minor;
l5 = legend('Buoyancy perturbation','Vertical velocity','Horizontal vorticity $\zeta$','interpreter','latex','Position',[0.07 0.2265 0.2063 0.0458]);
% l5 = legend('Buoyancy perturbation','Vertical velocity','interpreter','latex','Position',[0.07 0.2265 0.2063 0.0458]);
xlabel('Time (tidal cycles)','interpreter','latex')
set(gca,'Fontsize',fontsize);
title('Normalized perturbations','interpreter','latex','Fontsize',fontsize+5);
