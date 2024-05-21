

clear;close all;

dz = 1;   
Hmax = 250;
Nr = round(Hmax/dz);
zz_wgrid = 0:dz:((Nr)*dz);
t1hour = 3600;
DIV = 1;
FigureIsVisible = true;
fontsize = 20;
load_colors;
XLIM = [0 120];
NTtide = 50;

%%% Calculate the growth rate
Shear_parm = ([0:0.1:2.0])*1e-3;
lambda_max_all = [
10000
1260
2000
2000
1580
2000
2000
2000
1000
1000
2510
2510
2510
2000
2000
1000
710
560
560
560
560]';

for s = 2:21
    Shear = Shear_parm(s)
    lambda_max = lambda_max_all(s);
    expdir =  ['/Users/ysi/MITgcm_shear_convec/instability/exps_linear/lambda' ...
        num2str(lambda_max) '/topo0_H250_N0.001_S' num2str(Shear) ...
        '_lambda' num2str(lambda_max) '/'];

    % load([expdir 'output_2.mat'],'uuu','www','re_buoy','tt')
    load([expdir 'output_2.mat'])

    Nt = length(uuu);
    plot_tidx = 1:50:Nt;

    fit_span = round(Nt/NTtide*4)+1:round(Nt/NTtide*10);
    
    
    TKE = 0.5*(uuu.^2+0.5*(www(:,1:Nr)+www(:,2:Nr+1)).^2);
    TPE = 0;
    KE_PE = TKE+TPE;
    
    KE_PE_zavg = mean(KE_PE,2)';
    xxplot = tt/t1hour;
    yyplot = log(KE_PE_zavg)/2;
    [pKE,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
    GrowthRate_KE(s) = pKE(1);
    pKE(1);
    [y_fit,delta_fit] = polyval(pKE,xxplot,S);
    
    b2 = mean(re_buoy.^2,2)';
    yyplot_b2 = log(b2)/2;
    [pb2,S_b2] = polyfit(xxplot(fit_span),yyplot_b2(fit_span),1); 
    GrowthRate_b2(s) = pb2(1);
    [y_fit_b2,delta_fit_b2] = polyval(pb2,xxplot,S_b2);
    
    
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
xlim(XLIM)

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
xlim(XLIM)

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
xlim(XLIM)



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
xlim(XLIM)


subplot(3,2,6)
pcolor(tt(plot_tidx)/t1hour,zz_wgrid,www(plot_tidx,:)');shading flat;colorbar;
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');xlabel('Time (hours)')
title('Vertical velocity w','Fontsize',fontsize+3);
set(gca,'color',gray);
% if(kx~=0)
    caxis([-1 1]*max(max((abs(www))))/DIV)
% end
xlim(XLIM)

saveas(h,[expdir 'fig5.png'])

end

shear_Floquet = Shear_parm;
save('GrowthRate_exps_linear.mat','shear_Floquet','GrowthRate_b2','GrowthRate_KE')

figure(2)
plot(shear_Floquet,GrowthRate_b2)
hold on;
plot(shear_Floquet,GrowthRate_KE)
