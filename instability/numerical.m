



%%% Background tidal velocity
Atide = Shear*zz;
Atide_wgrid = Shear*zz_wgrid;
Utide =cos(tt*omega)'.*Atide;
% Utide = repmat(cos(tt*omega)',[1 length(Atide)])...
%     .*repmat(Atide,[length(tt) 1])/U0;

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

close all;

%%%% Integrate non-dimensionalized b and zeta with time
%%%% 1st-order and 2nd-order centered difference
%%%% Fourth-order Runge-Kutta or Third-order Adams-Bashforth
%%%% or Euler forward scheme for time advancement

zspan = [0 Hmax];

%%%%%%%%%%%% B.C.-1 %%%%%%%%%%%%
zeta(1,1) = 0; zeta(1,Nr+1) = 0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for o=1:Nt-1

    if(rem(o,round(Nt/10))==0)
        Progress = o/Nt
    end

    t0 = tt(o);
    b0 = buoy(o,:);
    z0 = zeta(o,:);
    
    loop;
    psi(o,:) = p0;

    if(useRK4)
        RK4;
    elseif(useAB3)
        %%% Third-order Adams-Bashforth method %%%
        if (o <= 2)
            %%% Use RK4 for the first 2 time steps
            RK4;
            %%% Use Euler-forward for the first 2 time steps
            % buoy(o+1,:) = buoy(o,:) + dbdt(o,:)*dt;
            % zeta(o+1,:) = zeta(o,:) + dzetadt(o,:)*dt;
        else
            buoy(o+1,:) = buoy(o,:) + dt*( (23/12)*dbdt(o,:)    - (16/12)*dbdt(o-1,:)    + (5/12)*dbdt(o-2,:) );
            zeta(o+1,:) = zeta(o,:) + dt*( (23/12)*dzetadt(o,:) - (16/12)*dzetadt(o-1,:) + (5/12)*dzetadt(o-2,:) );
        end
    else
        %%% Euler forward %%%
        buoy(o+1,:) = dbdt(o,:)*dt;
        zeta(o+1,:) = dzetadt(o,:)*dt;
    end

    %%%%%%%%%%%% B.C.-2 %%%%%%%%%%%%
    %%% No-stress (free-slip)
    zeta(o+1,1) = 0; zeta(o+1,Nr+1) = 0; 

    %%% No total stress at the ocean bottom
    % zeta(o,1) = cos(t0*omega);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

re_psi = real(psi);   re_psi(re_psi==0)=NaN;
re_zeta = real(zeta); re_zeta(re_zeta==0)=NaN;
re_buoy = real(buoy); re_buoy(re_buoy==0)=NaN;
re_dbdz = real(dbdz);

re_bq1 = real(bq1);
re_bq2 = real(bq2);
re_bq3 = real(bq3);
re_bq4 = real(bq4);
re_bq5 = real(bq5);

re_zq1 = real(zq1);
re_zq2 = real(zq2);
re_zq3 = real(zq3);
re_zq4 = real(zq4);


dbuoydz = zeros(Nt,Nr);
for m = 2:Nr-1
    dbuoydz(:,m) = (re_buoy(:,m+1)-re_buoy(:,m-1))/dz;
end


plot_tidx = 1:10:Nt;

load_colors;

DIV = 1;
uuu = -real((psi(:,2:Nr+1)-psi(:,1:Nr))/dz);
www = real(1i*kx*psi);

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

decompose;
close all

fit_span = round(Nt/NTtide*2)+1:Nt-1;

% clear TKE TPE KE_PE KE_PE_zavg TKE1 TKE2 p S 
TKE = 0.5*(uuu.^2+0.5*(www(:,1:Nr)+www(:,2:Nr+1)).^2);
TPE = 0;
KE_PE = TKE+TPE;

KE_PE_zavg = mean(KE_PE,2)';
xxplot = tt/t1hour;
yyplot = log(KE_PE_zavg)/2;
[pKE,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
GrowthRate_KE(Nexp_lambda,Nexp_shear) = pKE(1);
pKE(1);
[y_fit,delta_fit] = polyval(pKE,xxplot,S);

b2 = mean(re_buoy.^2,2)';
yyplot_b2 = log(b2)/2;
[pb2,S_b2] = polyfit(xxplot(fit_span),yyplot_b2(fit_span),1); 
GrowthRate(Nexp_lambda,Nexp_shear) = pb2(1);
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


close all;
% save(outputname)

clear b0 b_wgrid b0_wgrid b_2 b_3 b_4 p0 p0_ugrid psi psi0 sol1 solinit ...
    z0 z_2 z_3 z_4 zeta dbdz ...
    d2bdz2 d2psidz2 d2zetadz2  dpsidz dUtidedz dzetadt ...
    bq1 bq2 bq3 bq4 bq5 zq1 zq2 zq3 zq4 ...
    k_1b k_1z k_2b k_2z k_3b k_3z k_4b k_4z h ...
    buoy re_zq1 re_zq2 re_zq3 re_zq4 re_bq1 re_bq2 re_bq3 re_bq4 re_bq5 Utide ...
    uuu  re_dbdz re_d2bdz2 re_d2zetadz2 

outputname2 = [expdir 'output2.mat'];
save(outputname2)


