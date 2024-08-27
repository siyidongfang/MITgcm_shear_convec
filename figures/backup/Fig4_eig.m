
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
load('fig4/MIT_flat_noCori_tt.mat');
load('fig4/Floquet_km_flat_8e-5.mat')

% lam_x_real=lam_x_real/6;
% max_grow_rw = max(grow_rw,[],2);

for i=1:length(shear_MITgcm)
    [a(i) b(i)] = min(abs(shear_MITgcm(i)-shear_calc_Ri));
    Ri_gcm(i) = Ri_min(b(i));
end

for i=1:length(shear_floquet)
    [a(i) b(i)] = min(abs(shear_floquet(i)-shear_calc_Ri));
    Ri_floquet(i) = Ri_min(b(i));
end

for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km(i) = Ri_min(b(i));
end

% load('fig4/flat_km_kv2e-4.mat','max_grow','shear_all')
% for i=1:length(shear_all)
%     [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
%     Ri_km_diff(i) = Ri_min(b(i));
% end


clear shear_all
load('fig5/fig5_topo0_noDiff_eig.mat')
for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_eig(i) = Ri_min(b(i));
end



ax1 = subplot('position',[.065 .74 0.375 0.225]);
annotation('textbox',[0.028+0.02 0.996 0.15 0.01],'String','a','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
% pcolor(1./Ri_km,lam_x_real/1000,grow_rw');
pcolor(1./Ri_eig,m0_rw./kx_all,grow_all);

hold on;

xxxplot = 1./Ri_km;
xxxplot(xxxplot<1.97)=NaN;
yyyplot = N/omega.*sqrt(xxxplot-2);
plot(xxxplot,yyyplot,'LineWidth',2,'Color',black);
xxxplot = 1./Ri_km;
xxxplot(xxxplot<2.49)=NaN;
yyyplot = N/omega.*sqrt(xxxplot-2.5);
plot(xxxplot,yyyplot,'LineWidth',2,'Color',black);
xxxplot = 1./Ri_km;
xxxplot(xxxplot<2.99)=NaN;
yyyplot = N/omega.*sqrt(xxxplot-3);
plot(xxxplot,yyyplot,'LineWidth',2,'Color',black);
xxxplot = 1./Ri_km;
xxxplot(xxxplot<3.49)=NaN;
yyyplot = N/omega.*sqrt(xxxplot-3.5);
plot(xxxplot,yyyplot,'LineWidth',2,'Color',black);
xxxplot = 1./Ri_km;
yyyplot = N/omega.*sqrt(xxxplot+0.8);
plot(xxxplot,yyyplot,'LineWidth',2,'Color',black);
% xxxplot = 1./Ri_km;
% yyyplot = N/omega.*sqrt(xxxplot)+1;
% plot(xxxplot,yyyplot,'LineWidth',2,'Color',black);


ylim([0 1]*25)

% contour(1./Ri_km,1./rw_all,grow_rw',[0.02:0.02:0.5],'Color',black);


% set(gca, 'XScale', 'log')
shading interp;
set(gca,'Fontsize',fontsize);
colormap(WhiteBlueGreenYellowRed(0))
h1 = colorbar;
set(h1,'Position',[0.45 0.7456+0.01   0.008    0.2]);
set(get(h1,'Title'),'String',{'$\ \ \ \ (\mathrm{hour}^{-1})$'},'interpreter','latex','FontSize',fontsize);
clim([0 0.35]);
% xlabel('Minimum Richardson number ${R_i}_\mathrm{min}$','interpreter','latex');
xlabel('Inverse Richardson number ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
% ylabel('Horizontal wavelength (km)','interpreter','latex');
ylabel('Perturbation aspect ratio $m_0/k_0$','interpreter','latex');
title('Growth rate (flat bottom)','interpreter','latex','Fontsize',fontsize+5);
% ylim([0 33]/6)
xlim([0 4])
% xticks([0.1 0.2 0.3 0.5 0.75 1 2 3 4])


