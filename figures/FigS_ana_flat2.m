
clear;
addpath ../transition/
addpath ../analysis/
addpath ../analysis/colormaps/

load_colors;

fontsize = 18;

fg1 = figure(1);
clf;
set(gcf,'Color','w','Position', [100 153 1000 700])

tiledlay = tiledlayout(2,2);


load('../transition/freq_flat.mat')

nexttile
plot(R_all,om1,'LineWidth',2,'Color',green)
hold on;
plot(R_all,om2,'LineWidth',2,'Color',blue)
set(gca,'YScale', 'log');
grid on;grid minor;set(gca,'Fontsize',fontsize);
ylim([0.1 10])
xlabel('$R={R_i}_\mathrm{min}^{-1}$','Interpreter','latex');
legend('$\omega_1$ = harmonic mean of $\sigma(\hat t)\geq 1$',...
    '$\omega_2$ = harmonic mean of $\sigma(\hat t)<1$','Interpreter','latex','Fontsize',fontsize+2);
title('Idealized frequencies $\omega_1$ and $\omega_2$','Interpreter','latex','Fontsize',fontsize+4);
xlim([0 10])

nexttile
plot(R_all,alpha,'LineWidth',2,'Color',brown)
hold on;
plot(R_all,0.5*ones(1,length(R_all)),'--','LineWidth',2,'Color',yellow)
hold off;
grid on;grid minor;set(gca,'Fontsize',fontsize);
xlabel('$R={R_i}_\mathrm{min}^{-1}$','Interpreter','latex');
title('Portion of a period with $\sigma(\hat t)<1$','Interpreter','latex','Fontsize',fontsize+4);
yticks([0 0.5 1])
xlim([0 10])
ylim([0 1])

nexttile
plot(R_all,sqrt(A_harmonic),'LineWidth',2,'Color',brown)
hold on;
plot(R_all,ones(1,length(R_all)),'--','LineWidth',2,'Color',yellow)
grid on;grid minor;set(gca,'Fontsize',fontsize);
xlabel('$R={R_i}_\mathrm{min}^{-1}$','Interpreter','latex');
title('Harmonic mean of the frequency $\sigma(\hat t)$','Interpreter','latex','FontSize',fontsize+4);
legend('$\overline{\sigma}^\mathrm{hm}$','Unity','Interpreter','latex','FontSize',fontsize+4);
set(gca,'YScale', 'log');
ylim([0.1 10])
xlim([0 10])



load('../transition/grow_flat.mat')

load('../instability_km/exps_varyingN/N1e-3output','grow_rw','shear_all')
load('../figures/fig4/Ri_flat.mat')

max_grow_km = max(grow_rw,[],2);
for i=1:length(shear_all)
    [a(i) b(i)] = min(abs(shear_all(i)-shear_calc_Ri));
    Ri_km(i) = Ri_min(b(i));
end
grow_mzero = grow_rw(:,1);

nexttile
plot(R_all,grow,'LineWidth',2,'Color',black)
hold on;
plot(1./Ri_km,grow_mzero,'--','LineWidth',2)
xlim([0 10])
ylim([-0.01 0.4])
hold on;
grid on;grid minor;
xlabel('R = ${R_i}_\mathrm{min}^{-1}$','interpreter','latex');
set(gca,'Fontsize',fontsize);
legend('Analytical solution with $\omega_1$ and $\omega_2$','Inviscid theory ($m_0/k_0=0$)','Position',[0.5437 0.2961 0.3453 0.0750],'interpreter','latex','Fontsize',fontsize+2);
legend('boxoff')
title('Growth rate (1/hour)','Interpreter','latex','Fontsize',fontsize+4);


tiledlay.TileSpacing = 'compact';
tiledlay.Padding = 'compact';
AddLetters2Plots(fg1,'FontSize',fontsize+5,'FontWeight','normal')

print('-dpng','-r300',['fig_supp_new/figS_ana_flat2_matlab.png']);

