
clear; close all;
addpath ../transition/
addpath ../analysis/colormaps/

load_colors;

fontsize = 18;

fg1 = figure(1);
clf;
set(gcf,'Color','w','Position', [100 153 1200 800])

tiledlay = tiledlayout(2,2);

load('fig4/Ri_topo4_new.mat')
load('../instability_km/exps_new/topo4_nu0_large_shearoutput.mat')
for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km(i) = Ri_min(b(i));
end

r = cosd(topo)/sind(topo)-shear_all/omega;
for i=1:length(r)
    [a(i) b(i)] = min(abs(rw_all-1./r(i)));
    mr = rw_all(b(i));
    grow_mr(i) = grow_rw(i,b(i));
end

clear r_mostunstable
[max_grow r_idx]=max(grow_rw,[],2);
for i=1:length(shear_all)
    r_mostunstable(i) = 1./rw_all(r_idx(i));
end
% r_mostunstable(r_mostunstable>12)=NaN;


nexttile
rw_idx = 1:Nrw;
% grow_rw(:,1./rw_all>13)=NaN;
pcolor(1./Ri_km,1./rw_all,grow_rw');
hold on;
l11 = scatter(1./Ri_km,r_mostunstable,7,brown,'filled','LineWidth',1);
l12 = plot(1./Ri_km,r,'Color','k','LineWidth',2);
plot(1./Ri_km,zeros(1,length(Ri_km)),'k--','LineWidth',0.5);
set(gca,'Fontsize',fontsize);
shading interp;
colormap(WhiteBlueGreenYellowRed(0))
h2 = colorbar;
% set(h2,'Position',[0.45 0.4056+0.01   0.008    0.2]);
% set(get(h2,'Title'),'String',{'$\ \ \ \ (\mathrm{hour}^{-1})$'},'interpreter','latex','FontSize',fontsize);
clim([0 0.4]);
xlabel('${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
ylabel('Aspect ratio $r=m_0/k_0$','interpreter','latex');
title('Growth rate (1/hour)','interpreter','latex','Fontsize',fontsize+5);
ylim([-10 14])
xlim([0 8])
set(gca,'TickDir','out');
grid on;grid minor;
box on;
legend([l11 l12],'Diagnosed most unstable mode',...
    ['Predicted most unstable mode for ${R_i}_\mathrm{min}^{-1}>2$:' newline '$r=\cot\theta-\hat N \sqrt{R} = \cot\theta-\Lambda/\omega$'] , ...
    'interpreter','latex','Fontsize',fontsize+2,'Position',[0.0764 0.5813 0.2667 0.1220]);
legend('boxoff')




max_grow_rw = max(grow_rw,[],2);


nexttile
l21 = plot(1./Ri_km,max_grow_rw,'-','LineWidth',2,'Color',black);
grid on;grid minor;
hold on;
l22 = plot(1./Ri_km,grow_mr,'--','LineWidth',2,'Color',orange);
ylabel('(hour$^{-1}$)','interpreter','latex');
xlabel('${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
legend([l21 l22],'Growth rate of diagnosed most unstable modes','Growth rate along $r=\cot\theta-\Lambda/\omega$','interpreter','latex','Fontsize',fontsize+2,'Position', [0.5880 0.8830 0.3787 0.0666]);
legend('boxoff')

set(gca,'Fontsize',fontsize);
xlim([0 8])
ylim([-1e-3 0.35])
title('Growth rate','interpreter','latex','Fontsize',fontsize+5);
box on;





load('grow_topo_3om.mat')
grow_ana = grow;
shear_ana = sqrt(R_all.*N^2);


% load('/Users/ysi/MITgcm_shear_convec/instability_km/exps_new/topo4_nu0_output.mat','grow_rw','shear_all')
load('../instability_km/exps_new/topo4_nu0_output.mat','rw_all','grow_rw','shear_all')
load('../figures/fig4/Ri_topo4_new.mat')

for i=1:length(shear_ana)
    [a(i) b(i)] = min(abs(shear_ana(i)-shear_calc_Ri));
    Ri_ana(i) = Ri_min(b(i));
end

max_grow_km = max(grow_rw,[],2);
% for i=1:length(shear_all)
%     [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
%     Ri_km(i) = Ri_min(b(i));
% 
%     r=cosd(topo)/sind(topo)-N_hat*sqrt(1/Ri_km(i))
% 
%     [a(i) b(i)] = min(abs(r-1./rw_all));
%     grow_mr(i) = grow_rw(i,b(i));
% end



nexttile
plot(1./Ri_ana,grow_ana,'LineWidth',2,'Color',black)
hold on;
plot(1./Ri_km,grow_mr,'--','LineWidth',2)
xlim([0 8])
ylim([-1e-3 0.35])
xlabel('${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
legend('Analytical solution with $\omega_1$, $\omega_2$, and $\omega_3$','Inviscid theory: along $r=\cot\theta-\Lambda/\omega$','interpreter','latex','Fontsize',fontsize+3,...
    'Position',[0.0936 0.3793 0.3339 0.0686]);
legend('boxoff')

grid on;grid minor;set(gca,'Fontsize',fontsize);
title('Growth rate','Interpreter','latex','Fontsize',fontsize+4);
ylabel('(hour$^{-1}$)','interpreter','latex');


load('../transition/real_component.mat')
nexttile
plot(tt_hat/2/pi,sigma1_real,'LineWidth',2,'Color',blue)
grid on;grid minor;set(gca,'Fontsize',fontsize);
title('Real component of $\mu(\hat t)$ (Independent of ${R_i}_\mathrm{min}^{-1}$)','Interpreter','latex','Fontsize',fontsize+4);
xlabel('Time, $\hat t/(2\pi)$ (Tidal phase)','Interpreter','latex');



tiledlay.TileSpacing = 'compact';
tiledlay.Padding = 'compact';
AddLetters2Plots(fg1,'FontSize',fontsize+5,'FontWeight','normal')

print('-dpng','-r300',['fig_supp_new/figS_ana_topo2_matlab.png']);


