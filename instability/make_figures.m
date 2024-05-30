
FigureIsVisible = true;
fontsize = 20;
plot_tidx = 1:10:Nt;        
DIV = 1;
load_colors;

h=figure(1);
        set(h,'Visible', FigureIsVisible);clf;
        plot(Atide,zz,'LineWidth',2);
        title('Atide (m/s)');set(gca,'Fontsize',fontsize)
        grid on;grid minor;
        saveas(h,[expdir 'fig1.png'])
        
h=figure(2);
        set(h,'Visible', FigureIsVisible);clf;
        plot(diff(Atide)./diff(zz),0.5*(zz(1:end-1)+zz(2:end)),'LineWidth',2);
        title('Shear (1/s)');set(gca,'Fontsize',fontsize)
        grid on;grid minor;
        saveas(h,[expdir 'fig2.png'])
        
h=figure(3);
        set(h,'Visible', FigureIsVisible);clf;
        pcolor(tt/3600,zz,Utide');shading flat;colormap redblue; colorbar;
        title('Atide (m/s)');xlabel('Time (hours)');set(gca,'Fontsize',fontsize)
        saveas(h,[expdir 'fig3.png'])
        
        
h=figure(5);
        set(h,'color','w','Visible', FigureIsVisible,'Position',[67 346 1015 619]);
        clf;
        subplot(3,2,1)
        pcolor(tt(plot_tidx)/t1hour,zz_wgrid,re_psi(plot_tidx,:)');shading flat;colorbar;
        colormap(redblue)
        set(gca,'Fontsize',fontsize);
        ylabel('HAB (m)');
        title('Streamfunction \psi','Fontsize',fontsize+3);
        set(gca,'color',gray);
        aaa = max(max(abs(re_psi))/DIV);
        if(~isnan(aaa))
        caxis([-1 1]*aaa)
        end
        
        subplot(3,2,2)
        pcolor(tt(plot_tidx)/t1hour,zz_wgrid,re_zeta(plot_tidx,:)');shading flat;colorbar;
        set(gca,'Fontsize',fontsize);
        ylabel('HAB (m)');
        title('Horizontal vorticity perturbation \zeta','Fontsize',fontsize+3);
        set(gca,'color',gray);
        aaa = max(max(abs(re_zeta))/DIV);
        if(~isnan(aaa))
        caxis([-1 1]*aaa)
        end
        
        subplot(3,2,3)
        pcolor(tt(plot_tidx)/t1hour,zz,re_buoy(plot_tidx,:)');shading flat;colorbar;
        set(gca,'Fontsize',fontsize);
        ylabel('HAB (m)');
        title('Buoyancy perturbation b^\prime','Fontsize',fontsize+3);
        set(gca,'color',gray);
        aaa = max(max(abs(re_buoy))/DIV);
        if(~isnan(aaa))
        caxis([-1 1]*aaa)
        end
        
        
        subplot(3,2,5)
        pcolor(tt(plot_tidx)/t1hour,zz,uuu(plot_tidx,:)');shading flat;colorbar;
        colormap(redblue)
        set(gca,'Fontsize',fontsize);
        ylabel('HAB (m)');xlabel('Time (hours)')
        title('Horizontal velocity perturbation u^\prime','Fontsize',fontsize+3);
        set(gca,'color',gray);
        aaa = max(max((abs(uuu))))/DIV;
        if(~isnan(aaa)&& aaa~=0)
        caxis([-1 1]*aaa)
        end
        
        subplot(3,2,6)
        pcolor(tt(plot_tidx)/t1hour,zz_wgrid,www(plot_tidx,:)');shading flat;colorbar;
        set(gca,'Fontsize',fontsize);
        ylabel('HAB (m)');xlabel('Time (hours)')
        title('Vertical velocity w','Fontsize',fontsize+3);
        set(gca,'color',gray);
        if(kx~=0)
            caxis([-1 1]*max(max((abs(www))))/DIV)
        end
        
        saveas(h,[expdir 'fig5.png'])


xxplot = tt/t1hour;
yyplot = log(KE_zavg)/2;
yyplot_b2 = log(b2_zavg)/2;

h=figure(8);
        clf;
        set(h,'color','w','Visible', FigureIsVisible,'Position',[85 222 979 420]);
        plot(xxplot/12,yyplot,'LineWidth',2)
        hold on
        plot(xxplot/12,yyplot_b2,'LineWidth',2)
        plot(xxplot(fit_span)/12,y_fit(fit_span),':','LineWidth',2)
        plot(xxplot(fit_span)/12,y_fit_b2(fit_span),':','LineWidth',2)
        grid on;grid minor;
        set(gca,'Fontsize',20);
        % ylim([pKE(2)-3 pKE(2)+pKE(1)*max(xxplot)+2])
        xlabel('$t$ (tidal cycle)','Interpreter','Latex')
        ylabel('$\ln(e)/2$','Interpreter','Latex')
        hold off;axis tight
        legend('TKE','(b^\prime)^2','Position',[0.8141 0.1988 0.0684 0.1393])
        saveas(h,[expdir 'KE.png'])



pposition =  [153 197 1153 319];
blegend = [0.3120 0.1819 0.1425 0.3858];
zlegend =[0.7692 0.1934 0.1168 0.3031];

% blegend1 = [0.1351 0.1662 0.1425 0.3858];
blegend1 = [0.0761 0.0095 0.1425 0.3858];
zlegend1 = [0.5741 0.1558 0.1168 0.3031];

legend_b = {'$-Ub^\prime_x$',...
    '$-w^\prime \tilde N^2 \cos\theta$',...
    '$-u^\prime \tilde N^2 \sin\theta$',...
    '$-w^\prime B_z$',...
    'Diffusion'};

legend_zeta = {'$-U \zeta^\prime_x$',...
    '$b^\prime_x\cos\theta$',...
    '$-b^\prime_z\sin\theta$',...
    'Diffusion'};

lw = 2;

h=figure(6);
clf;
set(h,'color','w','Position',pposition,'Visible', FigureIsVisible);
subplot(1,2,1)
plot(tt/t1hour,bq1_int,'LineWidth',lw);
hold on;
plot(tt/t1hour,bq2_int,'LineWidth',lw);
plot(tt/t1hour,bq3_int,'LineWidth',lw);
plot(tt/t1hour,bq4_int,'LineWidth',lw);
plot(tt/t1hour,bq5_int,'LineWidth',lw,'Color',gray);
set(gca,'Fontsize',fontsize);
xlabel('Time (hours)')
legend(legend_b,'Interpreter','Latex','Fontsize',fontsize+3,'Position',blegend1);
grid on;grid minor;
title('Buoyancy budget')

subplot(1,2,2)
plot(tt/t1hour,(zq1_int),'LineWidth',lw);
hold on;
plot(tt/t1hour,(zq2_int),'LineWidth',lw);
plot(tt/t1hour,(zq3_int),'LineWidth',lw);
plot(tt/t1hour,(zq4_int),'LineWidth',lw,'Color',gray);
set(gca,'Fontsize',fontsize);xlabel('Time (hours)')
legend(legend_zeta,'Interpreter','Latex','Fontsize',fontsize+3,'Position',zlegend1);
grid on;grid minor;
title('Horizontal vorticity budget')

saveas(h,[expdir 'fig6.png'])


h=figure(7);
clf;
set(h,'color','w','Position',pposition,'Visible', FigureIsVisible);
subplot(1,2,1);
semilogy(tt/t1hour,abs(bq1_int),'LineWidth',lw);
hold on;
semilogy(tt/t1hour,abs(bq2_int),'LineWidth',lw);
semilogy(tt/t1hour,abs(bq3_int),'LineWidth',lw);
semilogy(tt/t1hour,abs(bq4_int),'LineWidth',lw);
semilogy(tt/t1hour,abs(bq5_int),'LineWidth',lw,'Color',gray);
set(gca,'Fontsize',fontsize);xlabel('Time (hours)')
legend(legend_b,'Interpreter','Latex','Fontsize',fontsize+3,'Position',blegend);
grid on;grid minor;
title('Buoyancy budget (absolute value, log axis)')

subplot(1,2,2)
semilogy(tt/t1hour,abs(zq1_int),'LineWidth',lw);
hold on;
semilogy(tt/t1hour,abs(zq2_int),'LineWidth',lw);
semilogy(tt/t1hour,abs(zq3_int),'LineWidth',lw);
semilogy(tt/t1hour,abs(zq4_int),'LineWidth',lw,'Color',gray);
set(gca,'Fontsize',fontsize);xlabel('Time (hours)')
legend(legend_zeta,'Interpreter','Latex','Fontsize',fontsize+3,'Position',zlegend);
grid on;grid minor;
title('Horizontal vorticity budget (absolute value, log axis)')

saveas(h,[expdir 'fig7.png'])

