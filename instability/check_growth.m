

load('exps_tanh_BottomCenter_dz2/lambda10000/topo0_H500_N0.001_S0.0005_lambda10000/output.mat')

fontsize  = 20;
fit_span = round(Nt/NTtide*2)+1:Nt-1;

% clear TKE TPE KE_PE KE_PE_zavg TKE1 TKE2 p S 
TKE = 0.5*(uuu.^2+0.5*(www(:,1:Nr)+www(:,2:Nr+1)).^2);
TPE = 0;
KE_PE = TKE+TPE;

KE_PE_zavg = mean(KE_PE,2)';
xxplot = tt/t1hour;
yyplot = log(KE_PE_zavg)/2;
[pKE,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
% GrowthRate_KE(Nexp_lambda,Nexp_shear) = pKE(1);
pKE(1);
[y_fit,delta_fit] = polyval(pKE,xxplot,S);

b2 = mean(re_buoy.^2,2)';
yyplot_b2 = log(b2)/2;
[pb2,S_b2] = polyfit(xxplot(fit_span),yyplot_b2(fit_span),1); 
% GrowthRate(Nexp_lambda,Nexp_shear) = pb2(1);
[y_fit_b2,delta_fit_b2] = polyval(pb2,xxplot,S_b2);


h=figure(8);
clf;
set(h,'color','w','Visible', true,'Position',[85 222 979 420]);
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
% saveas(h,[expdir 'KE.png'])



lw = 2;
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

load_colors;

h=figure(6);
clf;
set(h,'color','w','Position',pposition,'Visible', true);
subplot(1,2,1)
% plot(tt/t1hour,bq_all,'LineWidth',lw);
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
% semilogy(tt/t1hour,zq_all,'LineWidth',lw);
plot(tt/t1hour,(zq1_int),'LineWidth',lw);
hold on;
plot(tt/t1hour,(zq2_int),'LineWidth',lw);
plot(tt/t1hour,(zq3_int),'LineWidth',lw);
plot(tt/t1hour,(zq4_int),'LineWidth',lw,'Color',gray);
set(gca,'Fontsize',fontsize);xlabel('Time (hours)')
legend(legend_zeta,'Interpreter','Latex','Fontsize',fontsize+3,'Position',zlegend1);
grid on;grid minor;
title('Horizontal vorticity budget')

% saveas(h,[expdir 'fig6.png'])


h=figure(7);
clf;
set(h,'color','w','Position',pposition,'Visible', true);
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

% saveas(h,[expdir 'fig7.png'])

 zz_wgrid = 0:dz:((Nr)*dz);
plot_tidx = 1:10:Nt;
DIV = 1;
kx = 2*pi/lambda;
re_zeta = real(zeta);
h=figure(5);
set(h,'color','w','Visible', true,'Position',[67 346 1015 619]);
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

% saveas(h,[expdir 'fig5.png'])
