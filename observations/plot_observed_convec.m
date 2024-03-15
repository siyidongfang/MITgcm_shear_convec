%%

addpath /Users/ysi/MITgcm_shear_convec/analysis/
load_colors;

plot_tidx = 1:10:length(time_temp);


% figure(1);
% clf;set(gcf,'Color','w','Position',[114 662 1188*1.5 289]);
% subplot(1,3,1)
% pcolor(time_temp(plot_tidx),depth_temp,temp(plot_tidx,:)');shading flat;colorbar;
% hold on;
% contour(time_temp(plot_tidx),depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
% hold off;
% xlabel('Time (days)')
% ylabel('Depth (m)')
% set(gca,'Fontsize',fontsize);
% axis ij;
% clim([meanT-1 meanT+2])
% title('Temperature (^oC)')
% 
% subplot(1,3,2)
% pcolor(time_temp(plot_tidx),depth_temp,salt(plot_tidx,:)');shading flat;colorbar;
% hold on;
% contour(time_temp(plot_tidx),depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
% hold off;
% xlabel('Time (days)')
% ylabel('Depth (m)')
% set(gca,'Fontsize',fontsize);
% axis ij;
% clim([meanS-0.1 meanS+0.15])
% title('Estimated salinity (psu)')
% 
% subplot(1,3,3)
% pcolor(time_temp(plot_tidx),depth_n2,N2(plot_tidx,:)');shading flat;colorbar;
% hold on;
% contour(time_temp(plot_tidx),depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
% hold off;
% xlabel('Time (days)')
% ylabel('Depth (m)')
% set(gca,'Fontsize',fontsize);
% axis ij;
% title('N^2 (1/s^2)')
% clim([-1 1]/1e5)
% colormap(cmocean('balance'));


figure(2);
clf;set(gcf,'Color','w','Position', [669 91 718 416]);
subplot(1,2,1)
plot(T_tavg,depth_temp,'LineWidth',2);
axis ij;set(gca,'Fontsize',fontsize);
grid on;grid minor;
ylabel('Depth (m)')
title('Time mean temperature (^oC)')

subplot(1,2,2)
plot(N2_zavg,depth_n2,'LineWidth',2);
axis ij;set(gca,'Fontsize',fontsize);
grid on;grid minor;
ylabel('Depth (m)')
title('Time mean N^2 (1/s^2)')


%%


figure(4);
clf;set(gcf,'Color','w');
set(gcf,'Position', [56 352 1865 305]);
subplot(1,3,1)
pcolor(time_u,depth_u,uselect');shading flat;colorbar;colormap(redblue)
xlabel('Time (days)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([-0.3 0.3])
title('u (m/s)')

subplot(1,3,2)
pcolor(time_temp,depth_n2,uselect_n2grid');shading flat;colorbar;colormap(redblue)
xlabel('Time (days)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([-0.3 0.3])
title('u (m/s)')


% subplot(1,3,2)
% pcolor(time_u,depth_u,vselect');shading flat;colorbar;colormap(redblue)
% xlabel('Time (days)')
% ylabel('Depth (m)')
% set(gca,'Fontsize',fontsize);
% axis ij;
% clim([-0.3 0.3])
% title('v (m/s)')
% 
% subplot(1,3,3)
% pcolor(time_u,depth_u,wselect');shading flat;colorbar;colormap(redblue)
% xlabel('Time (days)')
% ylabel('Depth (m)')
% set(gca,'Fontsize',fontsize);
% axis ij;
% clim([-0.3 0.3]/6)
% title('w (m/s)')

%%


figure(9)
clf;set(gcf,'Color','w','Position', [56 352 1865 305]);
subplot(1,3,1)
pcolor(time_temp(plot_tidx),depth_n2,buoy1(plot_tidx,:)');shading flat;colorbar;
hold on;
contour(time_temp(plot_tidx),depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
hold off;
xlabel('Time (days)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
title('Buoyancy (m/s^2)')
colormap(WhiteBlueGreenYellowRed(0))
clim([-0.01 0.01]/4)

% subplot(1,3,2)
% pcolor(time_temp(plot_tidx),depth_n2,buoy2(plot_tidx,:)');shading flat;colorbar;
% hold on;
% contour(time_temp(plot_tidx),depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
% hold off;
% xlabel('Time (days)')
% ylabel('Depth (m)')
% set(gca,'Fontsize',fontsize);
% axis ij;
% title('Buoyancy (m/s^2)')
% colormap(WhiteBlueGreenYellowRed(0))
% clim([-0.01 0.01]/4)

% subplot(1,3,3)
% pcolor(time_temp(plot_tidx),depth_n2,buoy3(plot_tidx,:)');shading flat;colorbar;
% hold on;
% contour(time_temp(plot_tidx),depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
% hold off;
% xlabel('Time (days)')
% ylabel('Depth (m)')
% set(gca,'Fontsize',fontsize);
% axis ij;
% title('Buoyancy (m/s^2)')
% colormap(WhiteBlueGreenYellowRed(0))
% clim([-0.01 0.01]/4)




figure(10)
clf;set(gcf,'Color','w','Position', [56 352 1865 305]);
subplot(1,3,1)
pcolor(time_temp(plot_tidx),depth_reconst_n,n2_1(plot_tidx,:)');shading flat;colorbar;
hold on;
contour(time_temp(plot_tidx),depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
hold off;
xlabel('Time (days)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
title('N^2 (1/s^2)')
colormap(redblue)
clim([-1 1]/1e4/4)

% subplot(1,3,2)
% pcolor(time_temp(plot_tidx),depth_reconst_n,n2_2(plot_tidx,:)');shading flat;colorbar;
% hold on;
% contour(time_temp(plot_tidx),depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
% hold off;
% xlabel('Time (days)')
% ylabel('Depth (m)')
% set(gca,'Fontsize',fontsize);
% axis ij;
% title('N^2 (1/s^2)')
% clim([-1 1]/1e4/4)

% subplot(1,3,3)
% pcolor(time_temp(plot_tidx),depth_reconst_n,n2_3(plot_tidx,:)');shading flat;colorbar;
% hold on;
% contour(time_temp(plot_tidx),depth_temp,temp(plot_tidx,:)',meanT-2:0.5:meanT+2,'Color',black);
% hold off;
% xlabel('Time (days)')
% ylabel('Depth (m)')
% set(gca,'Fontsize',fontsize);
% axis ij;
% title('N^2 (1/s^2)')
% clim([-1 1]/1e4/4)

% figure(11)
% plot(N2(1,:),depth_n2); axis ij;
% hold on;
% plot(n20,depth_reconst_n); axis ij;


