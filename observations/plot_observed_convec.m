%%

addpath /Users/ysi/MITgcm_shear_convec/analysis/
load_colors;



% figure(1);
% clf;set(gcf,'Color','w','Position',[114 662 1188 289]);
% subplot(1,2,1)
% pcolor(time_temp,depth_temp,temp');shading flat;colorbar;
% hold on;
% contour(time_temp,depth_temp,temp',[meanT-2:0.5:meanT+2],'Color',black);
% hold off;
% xlabel('Time (hours)')
% ylabel('Depth (m)')
% set(gca,'Fontsize',fontsize);
% axis ij;
% clim([meanT-2 meanT+2])
% title('Temperature (^oC)')
% 
% subplot(1,2,2)
% pcolor(time_temp,depth_n2,N2');shading flat;colorbar;
% hold on;
% contour(time_temp,depth_temp,temp',meanT-2:0.5:meanT+2,'Color',black);
% hold off;
% xlabel('Time (hours)')
% ylabel('Depth (m)')
% set(gca,'Fontsize',fontsize);
% axis ij;
% title('N^2 (1/s^2)')
% clim([0 1.6]/1e5)
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
plot(N2avg,depth_n2,'LineWidth',2);
axis ij;set(gca,'Fontsize',fontsize);
grid on;grid minor;
ylabel('Depth (m)')
title('Time mean N^2 (1/s^2)')



figure(4);
clf;set(gcf,'Color','w');
set(gcf,'Position', [56 352 1865 305]);
subplot(1,3,1)
pcolor(time_uw,depth_uw,uselect');shading flat;colorbar;colormap(redblue)
xlabel('Time (hours)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([-0.3 0.3])
title('u (m/s)')

subplot(1,3,2)
pcolor(time_uw,depth_uw,vselect');shading flat;colorbar;colormap(redblue)
xlabel('Time (hours)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([-0.3 0.3])
title('v (m/s)')

subplot(1,3,3)
pcolor(time_uw,depth_uw,wselect');shading flat;colorbar;colormap(redblue)
xlabel('Time (hours)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([-0.3 0.3]/6)
title('w (m/s)')

%%

figure(6)
clf;set(gcf,'Color','w','Position',[114 662 1188 289]);
subplot(1,2,1)
pcolor(time_uw,depth_uw_stagger,temp_uw_stagger');shading flat;colorbar;
hold on;
contour(time_uw,depth_uw_stagger,temp_uw_stagger',meanT-2:0.5:meanT+2,'Color',black);
hold off;
xlabel('Time (hours)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
clim([meanT-2 meanT+2])
title('Temperature (^oC)')

subplot(1,2,2)
pcolor(time_uw,depth_uw,N2_uwgrid');shading flat;colorbar;
hold on;
contour(time_uw,depth_uw_stagger,temp_uw_stagger',meanT-2:0.5:meanT+2,'Color',black);
hold off;
xlabel('Time (hours)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
title('N^2 (1/s^2)')
clim([0 1.6]/1e5)
colormap(cmocean('balance'));


%%
% figure(7);
% clf;set(gcf,'Color','w','Position', [669 91 718 416]);
% subplot(1,2,1)
% plot(Ttavg_uw_stagger,depth_uw_stagger,'LineWidth',2);
% axis ij;set(gca,'Fontsize',fontsize);
% grid on;grid minor;
% ylabel('Depth (m)')
% title('2-day mean temperature (^oC)')
% 
% subplot(1,2,2)
% plot(N2avg_uwgrid,depth_uw,'LineWidth',2);
% axis ij;set(gca,'Fontsize',fontsize);
% grid on;grid minor;
% ylabel('Depth (m)')
% title('2-day mean N^2 (1/s^2)')

%%

time_uw = time_uw(ttselect);
temp_uw_stagger = temp_uw_stagger(ttselect,:);

figure(9)
clf;set(gcf,'Color','w','Position', [56 352 1865 305]);
subplot(1,3,1)
pcolor(time_uw,depth_uw,buoy1');shading flat;colorbar;
hold on;
contour(time_uw,depth_uw_stagger,temp_uw_stagger',meanT-2:0.5:meanT+2,'Color',black);
hold off;
xlabel('Time (hours)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
title('Buoyancy (m/s^2)')
colormap(WhiteBlueGreenYellowRed(0))
clim([0.005 0.015])

subplot(1,3,2)
pcolor(time_uw,depth_uw,buoy2');shading flat;colorbar;
hold on;
contour(time_uw,depth_uw_stagger,temp_uw_stagger',meanT-2:0.5:meanT+2,'Color',black);
hold off;
xlabel('Time (hours)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
title('Buoyancy (m/s^2)')
colormap(WhiteBlueGreenYellowRed(0))
clim([0.005 0.015])

subplot(1,3,3)
pcolor(time_uw,depth_uw,buoy3');shading flat;colorbar;
hold on;
contour(time_uw,depth_uw_stagger,temp_uw_stagger',meanT-2:0.5:meanT+2,'Color',black);
hold off;
xlabel('Time (hours)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
title('Buoyancy (m/s^2)')
colormap(WhiteBlueGreenYellowRed(0))
clim([0.005 0.015])


figure(10)
clf;set(gcf,'Color','w','Position', [56 352 1865 305]);
subplot(1,3,1)
pcolor(time_uw,depth_n2grid_recons,n2_1');shading flat;colorbar;
hold on;
contour(time_uw,depth_uw_stagger,temp_uw_stagger',meanT-2:0.5:meanT+2,'Color',black);
hold off;
xlabel('Time (hours)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
title('N^2 (1/s^2)')
colormap(redblue)
clim([-1 1]/1e4)

subplot(1,3,2)
pcolor(time_uw,depth_n2grid_recons,n2_2');shading flat;colorbar;
hold on;
contour(time_uw,depth_uw_stagger,temp_uw_stagger',meanT-2:0.5:meanT+2,'Color',black);
hold off;
xlabel('Time (hours)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
title('N^2 (1/s^2)')
clim([-1 1]/1e4)

subplot(1,3,3)
pcolor(time_uw,depth_n2grid_recons,n2_3');shading flat;colorbar;
hold on;
contour(time_uw,depth_uw_stagger,temp_uw_stagger',meanT-2:0.5:meanT+2,'Color',black);
hold off;
xlabel('Time (hours)')
ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);
axis ij;
title('N^2 (1/s^2)')
clim([-1 1]/1e4)




