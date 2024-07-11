


load('/Users/ysi/MITgcm_shear_convec/instability/exps_linear_dz0.5/lambda800/topo0_H250_N0.001_S0.0016_lambda800/output.mat')

fontsize = 15;
load_colors;

plot_tidx = Nt/25*14:Nt/25*17;  
DIV = 1.2;
xplot = tt(plot_tidx)/t1hour;
xplot = xplot-xplot(1);

figure(1)
clf;
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 800 600]);

%--- Initial temperature field, indicating the sponge layer
ax1 = subplot('position',[0.07 0.72 0.365 0.22]);
annotation('textbox',[0 0.99 0.15 0.01],'String','a','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(xplot,zz,re_buoy(plot_tidx,:)');shading flat;
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)','interpreter','latex');
title('Buoyancy perturbation $b^\prime$','Fontsize',fontsize+3,'interpreter','latex');
set(gca,'color',gray);
aaa = max(max(abs(re_buoy(plot_tidx,:)))/DIV);
clim([-1 1]*aaa)
xticks([0:6:36])
h1 = colorbar(ax1,'Position',[0.44 0.4+0.015 0.008 0.19]);
set(get(h1,'Title'),'String',{'$\ \ \ \ (\mathrm{1/s})$'},'interpreter','latex','FontSize',fontsize,'Position',[3.2000 128 0]);



ax2 = subplot('position',[0.57 0.72 0.365 0.22]);
annotation('textbox',[0.5 0.99 0.15 0.01],'String','b','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(xplot,zz_wgrid,re_zeta(plot_tidx,:)');shading flat;
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)','interpreter','latex');
title('Horizontal vorticity perturbation $\zeta$','Fontsize',fontsize+3,'interpreter','latex');
set(gca,'color',gray);
aaa = max(max(abs(re_zeta(plot_tidx,:)))/DIV);
clim([-1 1]*aaa*1.5)
xticks([0:6:36])
h2 = colorbar(ax2,'Position',[0.94 0.7200+0.015 0.008 0.19]);
set(get(h2,'Title'),'String',{'$\ \ \ \ (\mathrm{1/s})$'},'interpreter','latex','FontSize',fontsize,'Position',[3.2000 128 0]);


ax3 = subplot('position',[0.07 0.4 0.365 0.22]);
annotation('textbox',[0 0.665 0.15 0.01],'String','c','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(xplot,zz,uuu(plot_tidx,:)');shading flat;
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)','interpreter','latex');
title('Horizontal velocity perturbation $u^\prime$','Fontsize',fontsize+3,'interpreter','latex');
set(gca,'color',gray);
aaa = max(max((abs(uuu(plot_tidx,:)))))/DIV;
clim([-1 1]*aaa)
xticks([0:6:36])
h3 = colorbar(ax3,'Position',[0.94 0.4+0.015 0.008 0.19]);
set(get(h3,'Title'),'String',{'$\ \ \ \ (\mathrm{m/s})$'},'interpreter','latex','FontSize',fontsize,'Position',[3.2000 128 0]);




ax4 = subplot('position',[0.57 0.4 0.365 0.22]);
annotation('textbox',[0.5 0.665 0.15 0.01],'String','d','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(xplot,zz_wgrid,www(plot_tidx,:)');shading flat;
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)','interpreter','latex');
title('Vertical velocity perturbation $w^\prime$','Fontsize',fontsize+3,'interpreter','latex');
set(gca,'color',gray);
aaa = max(max((abs(www(plot_tidx,:)))))/DIV;
clim([-1 1]*aaa)
xticks([0:6:36])
h4 = colorbar(ax4,'Position',[0.44 0.08+0.015 0.008 0.19]);
set(get(h4,'Title'),'String',{'$\ \ \ \ (\mathrm{m/s})$'},'interpreter','latex','FontSize',fontsize,'Position',[3.2000 128 0]);



ax5 = subplot('position',[0.07 0.08 0.365 0.22]);
annotation('textbox',[0 0.345 0.15 0.01],'String','e','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(xplot,zz_wgrid,re_psi(plot_tidx,:)');shading flat;
colormap(cmocean('balanced'))
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)','interpreter','latex');
title('Streamfunction $\psi$','Fontsize',fontsize+3,'interpreter','latex');
set(gca,'color',gray);
aaa = max(max(abs(re_psi(plot_tidx,:)))/DIV);
clim([-1 1]*aaa)
xticks([0:6:36])
h5 = colorbar(ax5,'Position',[0.44 0.7200+0.015 0.008 0.19]);
set(get(h5,'Title'),'String',{'$\ \ \ \ (\mathrm{m^2/s})$'},'interpreter','latex','FontSize',fontsize,'Position',[3.2000 128 0]);
xlabel('Time (hours)','interpreter','latex');



max_ke = max(KE_zavg(plot_tidx));
max_b2 = max(b2_zavg(plot_tidx));
yyplot = log10(KE_zavg(plot_tidx)/max_ke);
yyplot_b2 = log10(b2_zavg(plot_tidx)/max_b2);

ax6 = subplot('position',[0.57 0.08 0.365 0.22]);
annotation('textbox',[0.5 0.345 0.15 0.01],'String','f','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
l2 = plot(xplot,yyplot_b2,'LineWidth',2,'Color','k');
hold on;
l1 = plot(xplot,yyplot,'LineWidth',2,'Color',green);

grid on;grid minor;
set(gca,'Fontsize',fontsize);
xlabel('Time (hours)','Interpreter','Latex')
ylabel('log(energy)','Interpreter','Latex')
hold off;axis tight
legend([l1 l2],'TKE','$(b^\prime)^2$','Interpreter','Latex','Position',[0.8331 0.0888 0.0921 0.0675],'Fontsize',fontsize)
xticks([0:6:36])
title('Normalized turbulent energy','Interpreter','Latex','Fontsize',fontsize+3)

print('-dpng','-r300',['fig_supp/figS_disc_timeseries.png']);