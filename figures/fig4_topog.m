
clear;close all;
addpath ../analysis/colormaps/
fontsize = 17;
load_colors;


figure(1)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 950 900]);


clear growth_MITgcm max_grow_rw Ri_gcm Ri_km shear_MITgcm shear_all shear_calc_Ri
%--- topo = 4 degrees
load('../instability_km/exps_new/topo4_nu0_output.mat')
load('fig4/Ri_topo4.mat')

lam_x_real = lam_x_real/6;
max_grow_rw = max(grow_rw,[],2);

% load('fig4/MIT_topo4_kv5e-6_tt.mat')

% for i=1:length(shear_MITgcm)
%     if(shear_MITgcm(i)>shear_convec)
%         shear_MITgcm(i)=NaN;
%         growth_MITgcm(i)=NaN;
%         Ri_gcm(i)=NaN;
%     else
%         [a(i) b(i)] = min(abs(shear_MITgcm(i)-shear_calc_Ri));
%         Ri_gcm(i) = Ri_min(b(i));
%     end
% end

for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km(i) = Ri_min(b(i));
end

r = cosd(topo)/sind(topo)-shear_all/omega;
[max_grow r_idx]=max(grow_rw,[],2);
for i=1:length(shear_all)
    r_mostunstable(i) = 1./rw_all(r_idx(i));
end

r_mostunstable(r_mostunstable>12)=NaN;

ax3 = subplot('position',[.06 .40 0.375 0.225]);
annotation('textbox',[0.028+0.02 0.65 0.15 0.01],'String','c','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
rw_idx = 1:Nrw;
pcolor(1./Ri_km,1./rw_all,grow_rw');
hold on;
plot(1./Ri_km,r,'Color','k','LineWidth',1.5);
scatter(1./Ri_km,r_mostunstable,8,brown,'filled','LineWidth',1);
hold off;
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
colormap(WhiteBlueGreenYellowRed(0))


