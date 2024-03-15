%%%
%%% calc_N2total.m
%%%
%%% Calculate the total stratification

load('/Users/ysi/MITgcm_shear_convec/instability/exps_RK4/lambda398/H300_topo4_Pt43200_N0.001_S0.0013_lambda398/output.mat')


dB0_dz = N^2 * cosd(topo) * ones(1,length(ttd));
dB_dz  = N^2 * (-sind(topo)*Shear/omega*sin(omega*ttd));

Utided = cos(omega*ttd)'.*Atide; % figure(3);pcolor(ttd/3600,zzd,Utided');shading flat;colormap redblue; colorbar;
re_buoyd = re_buoy*N^2*sind(topo)/omega;
% db_dx = 

bt1 = -www.*dB_dz';
bt2 = -www.*N^2*cosd(topo);
% bt3 = -Utided.*db_dx;

dtd = dt/omega;
dzd = dz*delta;

db_dz_w1 = diff(cumsum((bt1)*dt),1,2)/dzd;
db_dz_w2 = diff(cumsum((bt2)*dt),1,2)/dzd;


dbdz = diff(re_buoyd,1,2)/dzd;

zz_dbdz = 0.5*(zzd(1:end-1)+zzd(2:end));

pidx = round(Nt/30*4): round(Nt/30*12);
plot_time = ttd(pidx)/43200;
fontsize = 20;

figure(1)
pcolor(plot_time,zz_dbdz,dbdz(pidx,:)');shading flat;clim([-3 3]/1e6);colormap(redblue);colorbar


N2_total = dB0_dz'+dB_dz'+dbdz;

figure(2)
pcolor(plot_time,zz_dbdz,N2_total(pidx,:)');shading flat;clim([-3 3]/1e6);colormap(redblue);colorbar


db_dz_leveln = dbdz(:,30)';

figure(29)
pcolor(plot_time,zzd_wgrid,www(pidx,:)');shading flat;colormap redblue; colorbar;clim([-1 1]/10)


figure(30);clf
pcolor(plot_time,zzd,db_dz_w1(pidx,:)');shading flat;colormap redblue; colorbar;
clim([-1 1]/1e9)

figure(31);clf
pcolor(plot_time,zzd,db_dz_w2(pidx,:)');shading flat;colormap redblue; colorbar;
clim([-1 1]/1e9)

figure(4)
clf;set(gcf,'color','w');
hold on;grid on;grid minor;box on;
plot(plot_time,dB0_dz(pidx),'--','LineWidth',1.5);
plot(plot_time,dB_dz(pidx),'--','LineWidth',1.5);
plot(plot_time,db_dz_leveln(pidx),'--','LineWidth',1.5);
plot(plot_time,dB0_dz(pidx)+dB_dz(pidx)+db_dz_leveln(pidx),'LineWidth',2,'Color','k');
hold off;
set(gca,'FontSize',fontsize)
xlabel('Time (tidal cycles)');ylabel('N^2_{total} (s^{-2})')
xlim('tight')