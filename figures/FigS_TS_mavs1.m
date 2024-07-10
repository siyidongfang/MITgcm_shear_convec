clear;close all;
addpath /Users/ysi/Software/gsw_matlab_v3_06_11/thermodynamics_from_t/;
addpath /Users/ysi/Software/gsw_matlab_v3_06_11/library/;
addpath /Users/ysi/Software/gsw_matlab_v3_06_11/;
addpath ../analysis/colormaps/
addpath ../observations/
addpath ../observations/CTD/

fontsize = 17;
load_colors;

%%% make T-S diagram
load CTD.mat

tidx = 1:7;
lat_1to7 = mean(lat_all(tidx));
lon_1to7 = mean(lon_all(tidx));

pt_1to7 = pt_all(1:823,tidx);
CT_1to7 = CT_all(1:823,tidx);
psal_1to7 = psal_all(1:823,tidx);
SA_1to7 = SA_all(1:823,tidx);
pp_1to7 = P_all(1:823,tidx);

psal_1to7 = psal_1to7(:)';
pt_1to7 = pt_1to7(:)';
pp_1to7 = pp_1to7(:)';
p_ref = 0;
salt=psal_1to7;
pot_temp=pt_1to7;
depths=pp_1to7;
lat_ref=lat_1to7;
lon_ref=lon_1to7;

%%% Easier variable names
pt = pot_temp;
ss = salt;

pt_rest = pt_all(1:823,8:15);
ss_rest = psal_all(1:823,8:15);
pt_rest = pt_rest(:)';
ss_rest = ss_rest(:)';

%%% Get ranges of T/S
pt_max = max(max(pt)) + 0.5;
pt_min = min(min(pt)) - 0.5;
ss_max = max(max(ss)) + 0.05;
ss_min = min(min(ss)) - 0.05;

%%% Grid for contouring density
pt_step = (pt_max-pt_min)/100;
pt_grid = pt_min:pt_step:pt_max;
ss_step = (ss_max-ss_min)/100;
ss_grid = ss_min:ss_step:ss_max;
[PT_grid,SS_grid] = meshgrid(pt_grid,ss_grid);

%%% Use the GSW toolbox instead of densmdjwf
SA_grid = gsw_SA_from_SP(SS_grid,p_ref,lon_ref,lat_ref);  
CT_grid = gsw_CT_from_pt(SA_grid,PT_grid); 
pd = gsw_rho(SA_grid,CT_grid,p_ref) - 1000;


load('fig1/fig1.mat')
o1 = 0.078611022276845;
o2 = 34.638914089176580;
salt = o1*temp+o2;
plot_tidx = 1:length(time_temp);
meanS = mean(salt,'all','omitnan');

n2_1obs =  + cosd(topo)*(n2_1obs-meanN2*cosd(topo));
n2_1fit =  + cosd(topo)*(n2_1fit-meanN2*cosd(topo));


figure(1)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 1400 700]);

%--- T/S diagram
ax1 = subplot('position',[0.048 0.58 0.25 0.38]);
annotation('textbox',[0.028 0.993 0.15 0.01],'String','a','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
scatter(ss_rest,pt_rest,4,gray); 
hold on;
scatter(ss,pt,4,depths/1000); 
pot_dens_contours = [26.0:0.1:28.1];
[C,h] = contour(SS_grid,PT_grid,pd,pot_dens_contours,'EdgeColor','k');
% clabel(C,h);
clabel(C,h,'LabelSpacing',200);

min_depth = 1200-150;
min_depth = min_depth/1000;
max_depth = min([1400+150 depths(end)])/1000;
dp = 2/1000;
zidx = min_depth/dp:max_depth/dp;

[p,S] = polyfit(pt(zidx),ss(zidx),1); 
p1 = p(1);p2 = p(2);
plot(ss(zidx),ss(zidx)/p(1)-p(2)/p(1),'k--','LineWidth',2)
hold off;
xlabel('Practical salinity (psu)','interpreter','latex');
ylabel('Potential temperature ($^\circ$C)','interpreter','latex');
% clim([min(depths) max(depths)]);
clim([0 1700]/1000);
% axis([ss_min ss_max pt_min pt_max]);
axis([ss_min ss_max pt_min pt_max]);
box on;
set(gca,'FontSize',fontsize);
title('Temperature--salinity diagram','interpreter','latex','FontSize',fontsize+5);
% % colormap(flipdim(jet,2));
% mycolor=WhiteBlueGreenYellowRed(6);
% colormap(mycolor);
h1 = colorbar(ax1);
set(h1,'Position',[0.305  0.625 0.007 0.28]);
set(get(h1,'Title'),'String',{'$\ \ \ $Depth','$\ \ (\mathrm{km})$'},'Fontsize',fontsize,'interpreter','latex');



%--- salinity
ax2 = subplot('position',[0.38 0.58 0.25 0.38]);
annotation('textbox',[0.36 0.993 0.15 0.01],'String','b','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(time_temp(plot_tidx)*24,depth_temp,salt(plot_tidx,:)');
shading flat;
hold on;
contour(time_temp(plot_tidx)*24,depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
contour(time_temp(plot_tidx)*24,depth_n2,smooth_N2(plot_tidx,:)',[0 0],'Color',cyan,'LineWidth',0.75);
hold off;
xlabel('Time (hours)','interpreter','latex');
ylabel('Depth (m)','interpreter','latex');
set(gca,'Fontsize',fontsize);
axis ij;
% clim([meanS-0.1 meanS+0.15])
clim([min(min(salt))-0.02 max(max(salt))+0.02])
% clim([34.95 35.15]);
xlim([0 48])
colormap(cmocean('balance'))
title('Estimated salinity','Fontsize',fontsize+5,'interpreter','latex');
xticks([0:6:48]);
h2=colorbar(ax2);
set(h2,'Position',[0.635  0.625 0.007 0.28]);
set(get(h2,'Title'),'String','$\ \ \ \ (\mathrm{psu})$','Fontsize',fontsize,'interpreter','latex');
% freezeColors;



%--- N2
ax3 = subplot('position',[0.709 0.58 0.25 0.38]);
annotation('textbox',[0.69 0.993 0.15 0.01],'String','c','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(time_temp(plot_tidx)*24,depth_n2,smooth_N2(plot_tidx,:)');
hold on;
contour(time_temp(plot_tidx)*24,depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
contour(time_temp(plot_tidx)*24,depth_n2,smooth_N2(plot_tidx,:)',[0 0],'Color',cyan,'LineWidth',0.5);
hold off;
shading flat;
xlabel('Time (hours)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (m)','interpreter','latex','FontSize',fontsize);
set(gca,'Fontsize',fontsize);
axis ij;
clim([-1.75 1.75]/1e5)
colormap(cmocean('balance'))
title('Observed $\partial_{\tilde z} B$','Fontsize',fontsize+5,'interpreter','latex');
xticks([0:6:48]);
h3=colorbar(ax3);
ylim([min(depth_temp) max(depth_temp)]);xlim([0 48])
set(h3,'Position',[0.964  0.625 0.007 0.28]);
set(get(h3,'Title'),'String',{'$\ \ \ \ (1/\mathrm{s}^2)$',''},'Fontsize',fontsize,'interpreter','latex');
% freezeColors;



print('-dpng','-r300','fig_supp/figS_TS_mavs1_matlab.png');




