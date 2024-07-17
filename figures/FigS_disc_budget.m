
load_colors;
fontsize = 16;

load('/Users/ysi/MITgcm_shear_convec/instability/exps_topo4_linear_dz1/lambda700/topo4_H250_N0.001_S0.0014_lambda700/output.mat')

plot_tidx = Nt/10*3:Nt/10*7;  
xplot = tt(plot_tidx)/t1hour;
xplot = xplot-xplot(1);
XTICKS = [0:12:60];


legend_b = {'$-Ub^\prime_x$',...
    '$-w^\prime \tilde N^2 \cos\theta$',...
     '$-w^\prime B_z$',...
    '$-u^\prime \tilde N^2 \sin\theta$',...
    'Diffusion'};

legend_zeta = {'$-U \zeta^\prime_x$',...
    '$b^\prime_x\cos\theta$',...
    '$-b^\prime_z\sin\theta$',...
    'Diffusion'};


lw = 1.8;
FigureIsVisible = true;

figure(1)
clf;
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 800 600]);


mm1 = max(abs(bq1_int(plot_tidx)));
ax1 = subplot('position',[0.07 0.72-0.1 0.365 0.22+0.1]);
annotation('textbox',[0 0.99 0.15 0.01],'String','g','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
b5 = plot(xplot,bq5_int(plot_tidx)/mm1,'LineWidth',lw,'Color',gray);
hold on;
b1 = plot(xplot,bq1_int(plot_tidx)/mm1,'LineWidth',lw,'Color',blue);
b2 = plot(xplot,bq2_int(plot_tidx)/mm1,'LineWidth',lw,'Color',orange);
b4 = plot(xplot,bq4_int(plot_tidx)/mm1,'LineWidth',lw);
b3 = plot(xplot,bq3_int(plot_tidx)/mm1,'--','LineWidth',1,'Color',gold);
set(gca,'Fontsize',fontsize);
xlabel('Time (hours)','interpreter','latex');
ll1 = legend([b1 b2 b4 b3 b5],legend_b,'Interpreter','Latex','Position',[0.0745 0.6125 0.1457 0.1579]);
grid on;grid minor;
title('Buoyancy perturbation budget','Fontsize',fontsize+3,'interpreter','latex');
xticks(XTICKS)
xlim([0 48])


mm2 = max(abs(zq1_int(plot_tidx)));

ax2 = subplot('position',[0.57 0.72-0.1 0.365 0.22+0.1]);
annotation('textbox',[0.5 0.99 0.15 0.01],'String','h','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
z4 = plot(xplot,(zq4_int(plot_tidx))/mm2,'LineWidth',lw,'Color',gray);
hold on;
z1 = plot(xplot,(zq1_int(plot_tidx))/mm2,'LineWidth',lw,'Color',blue);
z2 = plot(xplot,(zq2_int(plot_tidx))/mm2,'LineWidth',lw,'Color',orange);
z3 = plot(xplot,(zq3_int(plot_tidx))/mm2,'--','LineWidth',1,'Color',gold);
set(gca,'Fontsize',fontsize);
xlabel('Time (hours)','interpreter','latex');
ll2 = legend([z1 z2 z3 z4],legend_zeta,'Interpreter','Latex','Position',[0.5774 0.6377 0.1216 0.1244]);
grid on;grid minor;
title('Horizontal vorticity perturbation budget','Fontsize',fontsize+3,'interpreter','latex');
xticks(XTICKS)
xlim([0 48])


print('-dpng','-r200',['fig_supp/FigS_disc_budget.png']);



