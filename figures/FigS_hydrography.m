clear;close all;
addpath ../analysis/colormaps/
fontsize = 16;
load_colors;


%%%------ Data for plotting bathymetry
% addpath ../observations/
% addpath ../observations/topography/
addpath ../observations/CTD/

load('fig_supp/FigS_hydrography.mat')

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



mvs1_lon = -11.861534;
mvs2_lon = -11.843511;

mvs1_lat = 54.198556;
mvs2_lat = 54.183718;

% [a b1] = min(abs(lonn-mvs1_lon))
[a b1] = min(abs(latn-mvs1_lat))
latn(b1);
lonn(b1);

% [a b2] = min(abs(lonn-mvs2_lon))
[a b2] = min(abs(latn-mvs2_lat))
latn(b2);
lonn(b2);


%--- Canyon depth and local topographic slope (in degrees and radians)
ax2 = subplot('position',[0.56 0.82 0.387 0.17]);
annotation('textbox',[0.495 0.99 0.15 0.01],'String','b','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
l1 = plot(along_canyonn,depthn,'-','LineWidth',3,'Color',black);
hold on;
s31 = scatter(along_canyonn(b1),depthn(b1),250,"+",'LineWidth',2,'MarkerEdgeColor',yellow); %%% MAVS 1&2
s32 = scatter(along_canyonn(b2),depthn(b2),250,"+",'LineWidth',2,'MarkerEdgeColor',yellow); %%% MAVS 1&2
l2 = plot(along_canyonn,-h_supercritical,'-','LineWidth',1.5,'Color',green);
hold off;
axis ij;grid on;grid minor;
set(gca,'FontSize',fontsize);
ylabel('Depth (m)','interpreter','latex');
xlim([0 max(along_canyonn)])
axis tight
leg2 = legend([l1 l2],'Canyon thalweg','Supercritical to the M2 tide',...
    'Position',[0.5672 0.9261 0.2406 0.0626],'interpreter','latex','Fontsize',fontsize+1);
legend boxoff;



ax3 = subplot('position',[0.56 0.61 0.387 0.17]);
annotation('textbox',[0.495 0.78 0.15 0.01],'String','c','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
yyaxis right
l31 = plot(along_midn,abs(s_topogn),'k-','LineWidth',1);
hold on;
s31 = scatter(along_midn(b1),abs(s_topogn(b1)),250,"+",'LineWidth',2,'MarkerEdgeColor',yellow); %%% MAVS 1&2
s32 = scatter(along_midn(b2),abs(s_topogn(b2)),250,"+",'LineWidth',2,'MarkerEdgeColor',yellow); %%% MAVS 1&2
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

leg3 = legend([l31 l32],'Topographic slope of the canyon thalweg','Critical slope for the M2 tide',...
    'Position',[0.5784 0.7216 0.3453 0.0626],'interpreter','latex','Fontsize',fontsize+1);
legend boxoff;


%%
%--- Potential temperature from CTD
ax4 = axes('position',[0.07 0.06 0.24 0.45]);
annotation('textbox',[0.02 0.53 0.15 0.01],'String','d','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');

h1 = plot(pt_all(:,1:7),P_all(:,1:7),'LineWidth',1,'Color',gray);
axis ij;grid on;grid minor;
hold on;
h2 = plot(pt_all(:,8:15),P_all(:,8:15),'LineWidth',1,'Color',gray);
h1_mean = plot(pt_mean15,pp15,'LineWidth',3,'Color','k');hold off;
xlabel('($^\circ$C)','interpreter','latex');
ylabel('Depth (m)','interpreter','latex');
set(gca,'Fontsize',fontsize)
title('Potential temperature','interpreter','latex','Fontsize',fontsize+4);
ylim([0 2700])
xlim([2 15])

%--- Salinity from CTD
ax5 = axes('position',[0.39 0.06 0.24 0.45]);
annotation('textbox',[0.34 0.53 0.15 0.01],'String','e','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
plot(psal_all(:,1:7),P_all(:,1:7),'LineWidth',1,'Color',gray);
axis ij;grid on;grid minor;
hold on;
plot(psal_all(:,8:15),P_all(:,8:15),'LineWidth',1,'Color',gray);
plot(psal_mean15,pp15,'LineWidth',3,'Color','k');hold off;
ylim([0 2700])
xlabel('(psu)','interpreter','latex');
ylabel('Depth (m)','interpreter','latex');
set(gca,'Fontsize',fontsize)
title('Salinity','interpreter','latex','Fontsize',fontsize+4);
xlim([34.9 35.45])

%--- N^2 from CTD
ax6 = axes('position',[0.7 0.06 0.24 0.45]);
annotation('textbox',[0.65 0.53 0.15 0.01],'String','f','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
plot(log10(N2_all(:,1:7)),Pmid_all(:,1:7),':','LineWidth',0.5,'Color',gray);
axis ij;grid on;grid minor;
hold on;
plot(log10(N2_all(:,8:15)),Pmid_all(:,8:15),':','LineWidth',0.5,'Color',gray);
N2_mean7(823)=NaN;
plot(log10(N2_mean15),Pmid_all(:,end),'LineWidth',1.5,'Color','k');
hold off;
ylim([0 2700])
xlim([-9 -3])
xlabel('$\log(\mathrm{s}^{-2})$','interpreter','latex');
ylabel('Depth (m)','interpreter','latex');
set(gca,'Fontsize',fontsize)
title('$\log(\partial_{\tilde z} b)$','interpreter','latex','Fontsize',fontsize+4);


% print('-dpng','-r300',['fig_supp/figS_hydrography_matlab.png']);