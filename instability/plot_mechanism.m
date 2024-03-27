%%%
%%% plot_mechanism.m
%%%
%%% Plot u, w, zeta, b of one tidal cycle

%%% Load the data
% load('/Users/ysi/MITgcm_shear_convec/instability/experiments/lambda100/H250_topo4_Pt43200_N0.001_S0.0015_lambda100/output.mat')

%%
%%% Select data for making plots
tidx = round(Nt/NTtide*7):round(Nt/NTtide*9);
b = re_buoyd(tidx,:);
u = uuu(tidx,:);
w = www(tidx,:);
z = re_zetad(tidx,:);
utide = Utide(tidx,:)*U0;
t = ttd(tidx)/3600;
t = t-t(1);



figure(1);clf;set(gcf,'Color','w','Position',[131 138 1452 524*1.6]);
subplot(3,2,1);
pcolor(t,zzd,b');
hold on;
contour(t,zzd,utide',[0.05:0.05:1],'color',darkgray);
contour(t,zzd,utide',[0 0],'color',darkgray,'LineWidth',1.5);
contour(t,zzd,utide',[-1:0.05:-0.05],'--','color',darkgray);
title('b^\prime (m/s^2)');
shading flat;colorbar;colormap(redblue);
clim([-1 1]*max(abs(b),[],'all')/2);
set(gca,'Fontsize',fontsize);
ylabel('HAB (m)');xlabel('Time (hour)');

subplot(3,2,2);
pcolor(t,zzd_wgrid,z');
hold on;
contour(t,zzd,utide',[0.05:0.05:1],'color',darkgray);
contour(t,zzd,utide',[0 0],'color',darkgray,'LineWidth',1.5);
contour(t,zzd,utide',[-1:0.05:-0.05],'--','color',darkgray);
title('\zeta^\prime (1/s)');shading flat;colorbar;colormap(redblue);
clim([-1 1]*max(abs(z),[],'all')/2);
set(gca,'Fontsize',fontsize);ylabel('HAB (m)');xlabel('Time (hour)');

subplot(3,2,3);
pcolor(t,zzd_wgrid,w');
hold on;
contour(t,zzd,utide',[0.05:0.05:1],'color',darkgray);
contour(t,zzd,utide',[0 0],'color',darkgray,'LineWidth',1.5);
contour(t,zzd,utide',[-1:0.05:-0.05],'--','color',darkgray);
title('w^\prime (m/s)');
shading flat;colorbar;colormap(redblue);
clim([-1 1]*max(abs(w),[],'all')/2);
set(gca,'Fontsize',fontsize);ylabel('HAB (m)');xlabel('Time (hour)');


subplot(3,2,4)
pcolor(t,zzd,u');
hold on;
contour(t,zzd,utide',[0.05:0.05:1],'color',darkgray);
contour(t,zzd,utide',[0 0],'color',darkgray,'LineWidth',1.5);
contour(t,zzd,utide',[-1:0.05:-0.05],'--','color',darkgray);
title('u^\prime (m/s)');
shading flat;colorbar;colormap(redblue);
clim([-1 1]*max(abs(u),[],'all')/2);
set(gca,'Fontsize',fontsize);ylabel('HAB (m)');xlabel('Time (hour)');



%%

TKE = 0.5*(u.^2+0.5*(w(:,1:Nr)+w(:,2:Nr+1)).^2);
TKE = TKE/(0.5*U0^2);
TPE = 0;
KE_PE = TKE+TPE;

KE_PE_zavg = mean(KE_PE,2)';
xxplot = t;
yyplot = log(KE_PE_zavg)/2;
b2 = mean(b.^2,2)';
yyplot_b2 = log(b2)/2;

yyplot_b2 = yyplot_b2 +8.5;

% h=figure(2);
% clf;
% set(gcf,'color','w','Position',[85 222 979 420]);
subplot(3,2,5)
plot(xxplot,yyplot,'LineWidth',2)
hold on
plot(xxplot,yyplot_b2,'LineWidth',2)
grid on;grid minor;
set(gca,'Fontsize',fontsize);
ylabel('$\ln(e)/2$','Interpreter','Latex')
xlabel('Time (hour)');
hold off;
legend('TKE','(b^\prime)^2','Position',[0.4105 0.1349 0.0475 0.0531])
axis tight
      



