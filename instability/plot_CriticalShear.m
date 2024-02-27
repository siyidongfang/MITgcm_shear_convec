
clear
fontsize = 20;
topo_parm = 0:0.01:20;
Ptide_parm = [1 2]*43200;

Ptide = Ptide_parm(1);
omega_m2 = 2*pi/Ptide;

Ptide = Ptide_parm(2);
omega_k1 = 2*pi/Ptide;

NT = length(topo_parm);
for nt = 1:NT
    topo = topo_parm(nt);    
    critical_convec_shear_m2(nt) = cotd(topo)*omega_m2;
    critical_convec_shear_k1(nt) = cotd(topo)*omega_k1;
end

topo_parm_plot = topo_parm/180*pi;
% topo_parm_plot = topo_parm;
figure(11)
clf;
set(gcf,'color','w','Position',[133 306 717 582]);
l1 = plot(topo_parm_plot,critical_convec_shear_m2,'LineWidth',2);
hold on;
sp1 =plot([0 topo_parm_plot],0.8e-3*ones(1,length(critical_convec_shear_m2)+1),...
    '--','LineWidth',2);
l2 = plot(topo_parm_plot,critical_convec_shear_k1,'LineWidth',2);

legend([l1 l2 sp1],'Semidiurnal tide','Diurnal tide','Observed shear')
grid on;grid minor;
set(gca,'Fontsize',fontsize);
ylabel('\Lambda_c (1/s)')
xlabel('Slope \theta')
title({'Critical shear of background stratification',...
    'crossing the zero line in one tidal period'})
ylim([0 5]*1e-3)

%%



clear;
fontsize = 20;

load('CriticalShear.mat')
% %%
% 
% figure(1)
% set(gcf,'color','w','Position',[133 306 717 582]);
% semilogy(Shear_parm,Ri_min,'LineWidth',2)
% hold on;plot(Shear_parm,0.25*ones(1,length(Shear_parm)),'-.','LineWidth',2)
% grid on;
% set(gca,'Fontsize',fontsize);
% ylabel('Ri_{min}')
% xlabel('Shear (1/s)')
% title('Minimum Ri of the background flow','Fontsize',fontsize+3);

load_colors;

%%
topo_parm_plot = topo_parm/180*pi;

figure(2)
clf;set(gcf,'color','w');
h=pcolor(log10(Nsquare_parm),topo_parm_plot,log10(criticalShear_M2)');
shading interp;
c=colorbar;
colormap(WhiteBlueGreenYellowRed(0));
hold on;
l1 = contour(log10(Nsquare_parm),topo_parm_plot,log10(criticalShear_M2)',[log10(0.8e-3) log10(0.8e-3)],'LineWidth',2,'Color',black);
l2 = plot(-6*ones(1,length(topo_parm_plot)),topo_parm_plot,'--','Color',purple);
xlabel('Background N^2 (1/s^2)')
ylabel('Slope \theta')
title({'Critical shear of background Ri \leq 1/4','\Lambda_s (1/s)'})
set(gca,'Fontsize',fontsize,...
    'XTick',[-9:1:-4],...
    'XTickLabel',{'10^{-9}','10^{-8}','10^{-7}','10^{-6}','10^{-5}','10^{-4}'});
c.Ticks=[-4 -3.5 -3 -2.5 -2];
c.TickLabels = {'10^{-4}','10^{-3.5}','10^{-3}','10^{-2.5}','10^{-2}'};

%%

zcolor=log10(Nsquare_parm);
cmap = colormap(WhiteBlueGreenYellowRed(8));
zmap = linspace( min(zcolor), max(zcolor), length(cmap));

figure(3)
clf;set(gcf,'color','w');
i=1;
color = interp1( zmap, cmap, zcolor(i));
loglog(Bu(i,:),criticalShear_M2(i,:),'LineWidth',1.5,'Color',color);
hold on;
for i=2:length(Nsquare_parm)
% for i=10
    color = interp1( zmap, cmap, zcolor(i));
    loglog(Bu(i,:),criticalShear_M2(i,:),'LineWidth',1.5,'Color',color);
end
loglog(0:0.005:70,0.8e-3*ones(1,length(0:0.005:70)),'LineWidth',1.5,'Color',black);
set(gca,'Fontsize',fontsize)
xlabel('Slope Burger number Bu = \theta N/f')
ylabel('\Lambda_s (1/s)')
title({'Critical shear of minimum background Ri \leq 1/4','for different background N^2'})
clim([min(zcolor), max(zcolor)])
c=colorbar;
c.Ticks=-9:1:-4;
c.TickLabels = {'10^{-9}','10^{-8}','10^{-7}','10^{-6}','10^{-5}','10^{-4}'};
grid on;


%%

figure(4)
clf;set(gcf,'color','w');
h=pcolor(log10(Nsquare_parm),topo_parm_plot,log10(Bu)');
shading interp;
c=colorbar;
colormap(WhiteBlueGreenYellowRed(0));
hold on;
% l1 = contour(log10(Nsquare_parm),topo_parm_plot,log10(Bu)',[log10(0.8e-3) log10(0.8e-3)],'LineWidth',2,'Color',black);
% l2 = plot(-6*ones(1,length(topo_parm_plot)),topo_parm_plot,'--','Color',purple);
xlabel('Background N^2 (1/s^2)')
ylabel('Slope \theta')
title('Slope Burger number Bu = \theta N/f')
set(gca,'Fontsize',fontsize,...
    'XTick',[-9:1:-4],...
    'XTickLabel',{'10^{-9}','10^{-8}','10^{-7}','10^{-6}','10^{-5}','10^{-4}'});
c.Ticks=[-2 -1 0 1];
c.TickLabels = {'0.01','0.1','1','10'};



%%


zcolor=topo_parm_plot;
cmap = colormap(WhiteBlueGreenYellowRed(8));
zmap = linspace( min(zcolor), max(zcolor), length(cmap));

figure(5)
clf;set(gcf,'color','w');
i=1;
color = interp1( zmap, cmap, zcolor(i));
loglog(Bu(:,i),criticalShear_M2(:,i),'LineWidth',1.5,'Color',color);
hold on;
for i=2:length(topo_parm)
% for i=10
    color = interp1( zmap, cmap, zcolor(i));
    loglog(Bu(:,i),criticalShear_M2(:,i),'LineWidth',1.5,'Color',color);
end
loglog(0:0.005:70,0.8e-3*ones(1,length(0:0.005:70)),'LineWidth',1.5,'Color',black);
set(gca,'Fontsize',fontsize)
xlabel('Slope Burger number Bu = \theta N/f')
ylabel('\Lambda_s (1/s)')
title({'Critical shear of minimum background Ri \leq 1/4','for different topographic slope \theta'})
clim([min(zcolor), max(zcolor)])
c=colorbar;
% c.Ticks=-9:1:-4;
% c.TickLabels = {'10^{-9}','10^{-8}','10^{-7}','10^{-6}','10^{-5}','10^{-4}'};
grid on;
