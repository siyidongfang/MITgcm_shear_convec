clear;close all;
addpath ../analysis/colormaps/
fontsize = 16;
load_colors;


%%%------ Data for plotting bathymetry
% addpath ../observations/
% addpath ../observations/topography/
% addpath ../observations/CTD/

load('fig_supp/FigS_obs_CTD.mat')

%%% Calculate the slope of internal wave characteristics r_iw (Lamb 2013, Annual Reviews, Eq. 3)
%%% Calculate the topographic slope s_topog
%%% Use the mean buoyancy frequency N 100m above the topography to calculate r_iw
%%% Compare r_iw with s_topog

omega_M2 = 2*pi/44712;
omega_idealized = 2*pi/43200;
latN = 54.2;
f0 = 4*pi/86164*sind(latN); 
load CTD.mat

omega_tides = omega_M2;
h = -depthn;
h_mid = 0.5*(h(1:end-1)+h(2:end));

Ny = length(depthn);
r_iw = NaN*zeros(1,Ny-1);


for i=1:Ny-1
  clear zidx_100 zidx1
  [a,zidx1] = min(abs(-h_mid(i)-p_mid15));
  zidx_100m = zidx1-50:zidx1;
  N_mean(i) = sqrt(mean(N2_mean15(zidx_100m)));
  r_iw(i) = sqrt((omega_tides^2-f0^2)/(N_mean(i)^2-omega_tides^2));
end


h_subcritical = NaN.*zeros(1,Ny);
h_supercritical = NaN.*zeros(1,Ny);

for i=1:Ny-1
  if(abs(s_topogn(i))<r_iw(i))
      h_subcritical(i) = h(i);
  else
      h_supercritical(i) = h(i);
  end
end


%%

figure(1)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 1000 800]);

%--- Canyon bathymetry, flat view
ax1 = subplot('position',[0.07 0.61 0.35 0.35]);
annotation('textbox',[0.04 0.99 0.15 0.01],'String','a','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
contour(lon,lat,z',[400:25:3000],'LineWidth',0.5);
hold on;
contour(lon,lat,z',[400:100:3000],'LineWidth',2,'ShowText','on');
hold off;
shading flat;
colormap(cmocean('rain'));
clim([800 2700])
set(gca,'FontSize',fontsize);
title('Canyon bathymetry (m)','FontSize',fontsize+4,'interpreter','latex')
xlabel('Longitude','interpreter','latex');
ylabel('Latitude','interpreter','latex');
xlim([-12.13 -11.8])
hold on;scatter(lonn,latn,60,".",'LineWidth',1,'MarkerEdgeColor',black);
s1 = scatter(lon_CTD(8:15),lat_CTD(8:15),200,"x",'LineWidth',4,'MarkerEdgeColor',orange);
scatter(lon_CTD(1:7),lat_CTD(1:7),200,"x",'LineWidth',4,'MarkerEdgeColor',orange);
s2 = scatter([-11.861534 -11.843511],[54.198556 54.183718],100,"+",'LineWidth',6,'MarkerEdgeColor',yellow); %%% MAVS 1&2
% scatter([-11-56.923/60 -11-52.268/60],[54+14.312/60 54+12.167/60],150,"^",'LineWidth',4); %%% MP 1&2
ylim([ 54.14 54.33])
xlim([-12.06 -11.8])
% scatter(lon_CTD(8:15),lat_CTD(8:15),150,"x",'LineWidth',4);
box on;grid on;grid minor;
hold off;
h1 = colorbar;
set(h1,'Position', [0.43 0.65 0.008 0.3], 'YDir', 'reverse');

leg1 = legend([s1 s2],'CTD','MAVS','Position',[0.0744 0.6221 0.0896 0.0566],...
    'interpreter','latex','Fontsize',fontsize+1);
legend boxoff;


%--- Canyon depth and local topographic slope (in degrees and radians)
ax2 = subplot('position',[0.56 0.82 0.387 0.17]);
annotation('textbox',[0.495 0.99 0.15 0.01],'String','b','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
l1 = plot(along_canyonn,depthn,'-','LineWidth',3,'Color',black);
hold on;
l2 = plot(along_canyonn,-h_supercritical,'-','LineWidth',1.5,'Color',green);
hold off;
axis ij;grid on;grid minor;
set(gca,'FontSize',fontsize);
ylabel('Depth (m)','interpreter','latex');
xlim([0 max(along_canyonn)])
axis tight
leg2 = legend([l1 l2],'Canyon thalweg','Supercritical to M2 tide',...
    'Position',[0.5672 0.9261 0.2406 0.0626],'interpreter','latex','Fontsize',fontsize+1);
legend boxoff;



ax3 = subplot('position',[0.56 0.61 0.387 0.17]);
annotation('textbox',[0.495 0.78 0.15 0.01],'String','c','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
yyaxis right
l31 = plot(along_midn,abs(s_topogn),'k-','LineWidth',1);
hold on;
l32 = plot(along_midn,abs(r_iw),'-.','LineWidth',2,'Color',green);
xlabel('Along-canyon distance (km)','interpreter','latex');
ylabel('(radian)','interpreter','latex');
ax3.YAxis(1).Color = black;
grid on;grid minor;
ylim([0 22/180*pi])

yyaxis left
plot(along_midn,abs(s_topogn)*180/pi,'k:','LineWidth',0.01);
ylabel('(degree)','interpreter','latex');
axis tight
ax3.YAxis(2).Color = black;
xlim([0 max(along_canyonn)])
ylim([0 22])
set(gca,'FontSize',fontsize);

leg3 = legend([l31 l32],'Topographic slope of the canyon thalweg','Critical slope for M2 tide',...
    'Position',[0.5784 0.7216 0.3453 0.0626],'interpreter','latex','Fontsize',fontsize+1);
legend boxoff;


%%
%--- Potential temperature from CTD
ax4 = axes('position',[0.07 0.05 0.24 0.45]);
annotation('textbox',[0.02 0.53 0.15 0.01],'String','d','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');




%--- Salinity from CTD

%--- N^2 from CTD


% print('-dpng','-r200',['fig_supp/figS_obs_CTD.png']);