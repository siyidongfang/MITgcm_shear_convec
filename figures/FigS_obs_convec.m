clear;close all;
addpath /Users/ysi/Software/gsw_matlab_v3_06_11/thermodynamics_from_t/;
addpath /Users/ysi/Software/gsw_matlab_v3_06_11/library/;
addpath /Users/ysi/Software/gsw_matlab_v3_06_11/;
addpath ../analysis/colormaps/
addpath ../observations/
addpath ../observations/CTD/

fontsize = 16;
load_colors;


load('fig1/fig1.mat')
plot_tidx = 1:1:length(time_temp);

zz_rec = -depth_n2 + depth_n2(end);
buoy_init = repmat(meanN2*cosd(topo)*zz_rec,[length(buoy1fit) 1]);
buoy_fit = buoy1fit + buoy_init;
buoy_obs = buoy1obs + buoy_init;

n2_1obs = meanN2 + cosd(topo)*(n2_1obs-meanN2*cosd(topo));
n2_1fit = meanN2 + cosd(topo)*(n2_1fit-meanN2*cosd(topo));

n2_1obs_convec =  cosd(topo)*(n2_1obs-meanN2*cosd(topo));
n2_1fit_convec =  cosd(topo)*(n2_1fit-meanN2*cosd(topo));



figure(1)
clf;   
set(gcf,'Color','w');
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0.03*scrsz(3) 0.3*scrsz(4) 1400 700]);


%%% Observed velocity
ax1 = subplot('position',[0.38 0.58-0.005 0.25 0.38]);
annotation('textbox',[0.36 0.993-0.005 0.15 0.01],'String','a','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(time_temp(plot_tidx)*24,depth_reconst_n,n2_1obs(plot_tidx,:)');
shading flat;
hold on;
contour(time_temp(plot_tidx)*24,depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
contour(time_temp(plot_tidx)*24,depth_reconst_n,n2_1obs(plot_tidx,:)',[0 0],'Color','c','LineWidth',1.5);
hold off;
xlabel('Time (hours)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (m)','interpreter','latex','FontSize',fontsize);
set(gca,'Fontsize',fontsize);
axis ij;
clim([-1.75 1.75]/1e5)
colormap(cmocean('balance'))
title('$\partial_{\tilde z} \mathcal B_\mathrm{rec}$ with $u_\mathrm{obs},\ \partial_{\tilde z}\mathcal B_\mathrm{rec}\big\vert_{t=0}={\tilde N}^2$','Fontsize',fontsize+5,'interpreter','latex');
xticks([0:6:48])
h1=colorbar(ax1);
ylim([min(depth_temp) max(depth_temp)]);xlim([0 48])
set(h1,'Position',[0.635  0.645-0.005 0.007 0.28]);
set(get(h1,'Title'),'String',{'$\ \ \ \ (1/\mathrm{s}^2)$',''},'Fontsize',fontsize,'interpreter','latex');

%%% Linear-fit velocity
ax2 = subplot('position',[0.709 0.58-0.005 0.25 0.38]);
annotation('textbox',[0.69 0.993-0.005 0.15 0.01],'String','b','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(time_temp(plot_tidx)*24,depth_reconst_n,n2_1fit(plot_tidx,:)');
shading flat;
hold on;
contour(time_temp(plot_tidx)*24,depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
contour(time_temp(plot_tidx)*24,depth_reconst_n,n2_1fit(plot_tidx,:)',[0 0],'Color','c','LineWidth',1.5);
hold off;
xlabel('Time (hours)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (m)','interpreter','latex','FontSize',fontsize);
set(gca,'Fontsize',fontsize);
axis ij;
clim([-1.75 1.75]/1e5)
colormap(cmocean('balance'))
title('$\partial_{\tilde z} \mathcal B_\mathrm{rec}$ with $u_\mathrm{fit},\ \partial_{\tilde z}\mathcal B_\mathrm{rec}\big\vert_{t=0}={\tilde N}^2$','Fontsize',fontsize+5,'interpreter','latex');
xticks([0:6:48])
ylim([min(depth_temp) max(depth_temp)]);xlim([0 48])
h2=colorbar(ax2);
set(h2,'Position',[0.964  0.645-0.005 0.007 0.28]);
set(get(h2,'Title'),'String',{'$\ \ \ \ (1/\mathrm{s}^2)$',''},'Fontsize',fontsize,'interpreter','latex');




%--- dbdz using observed u, without N2
ax5 = subplot('position',[0.38 0.07-0.005 0.25 0.38]);
annotation('textbox',[0.36 0.482-0.005 0.15 0.01],'String','c','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(time_temp(plot_tidx)*24,depth_reconst_n,n2_1obs_convec(plot_tidx,:)');
shading flat;
hold on;
contour(time_temp(plot_tidx)*24,depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
contour(time_temp(plot_tidx)*24,depth_reconst_n,n2_1obs_convec(plot_tidx,:)',[0 0],'Color','c','LineWidth',1.5);
hold off;
xlabel('Time (hours)','interpreter','latex','FontSize',fontsize);
ylabel('Depth (m)','interpreter','latex','FontSize',fontsize);
set(gca,'Fontsize',fontsize);
axis ij;
clim([-1.75 1.75]/1e5)
colormap(cmocean('balance'))
title('$\partial_{\tilde z} \mathcal B_\mathrm{rec}$ with $u_\mathrm{obs},\ \partial_{\tilde z}\mathcal B_\mathrm{rec}\big\vert_{t=0}=0$','Fontsize',fontsize+5,'interpreter','latex');
xticks([0:6:48])
h5=colorbar(ax5);
set(h5,'Position',[0.635 0.135-0.005 0.007 0.28]);
set(get(h5,'Title'),'String',{'$\ \ \ \ (1/\mathrm{s}^2)$',''},'Fontsize',fontsize,'interpreter','latex');
xlim([0 48])

%--- dbdz using linear-fit u, without N2
ax6 = subplot('position',[0.709 0.07-0.005 0.25 0.38]);
annotation('textbox',[0.69 0.482-0.005 0.15 0.01],'String','d','FontSize',fontsize+3,'fontweight','bold','LineStyle','None');
pcolor(time_temp(plot_tidx)*24,depth_reconst_n,n2_1fit_convec(plot_tidx,:)');shading flat;
hold on;
contour(time_temp(plot_tidx)*24,depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
contour(time_temp(plot_tidx)*24,depth_reconst_n,n2_1fit_convec(plot_tidx,:)',[0 0],'Color','c','LineWidth',1.5);
hold off;
xlabel('Time (hours)','interpreter','latex','FontSize',fontsize);
% ylabel('Depth (m)','interpreter','latex','FontSize',fontsize);
set(gca,'Fontsize',fontsize);
axis ij;
clim([-1.75 1.75]/1e5)
colormap(cmocean('balance'))
title('$\partial_{\tilde z} \mathcal B_\mathrm{rec}$ with $U_\mathrm{fit},\ \partial_{\tilde z}\mathcal B_\mathrm{rec}\big\vert_{t=0}=0$','Fontsize',fontsize+5,'interpreter','latex');
xticks([0:6:48])
h6=colorbar(ax6);
set(h6,'Position',[0.964 0.135-0.005 0.007 0.28]);
set(get(h6,'Title'),'String',{'$\ \ \ \ (1/\mathrm{s}^2)$',''},'Fontsize',fontsize,'interpreter','latex');
xlim([0 48])



print('-dpng','-r300','fig_supp/figS_obs_convec_matlab.png');
