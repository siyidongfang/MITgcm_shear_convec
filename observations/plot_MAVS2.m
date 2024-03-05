
addpath ../analysis/colormaps/
load MAVS2_Phase_Average.mat

fontsize = 20;

time = sin(phase);

figure(1)
clf;set(gcf,'color','w')
subplot(2,1,1)
pcolor(phase/pi*180,hab(:,1),temp)
shading flat;colorbar;colormap(redblue)
xlabel('Tidal phase (degrees)');ylabel('HAB (m)')
set(gca,'Fontsize',fontsize)
title('Temperature (^oC)','Fontsize',fontsize+2)

subplot(2,1,2)
pcolor(phase/pi*180,hab(:,1),N2)
shading flat;colorbar;colormap(redblue)
xlabel('Tidal phase (degrees)');ylabel('HAB (m)')
set(gca,'Fontsize',fontsize)
title('N^2 (1/s^2)','Fontsize',fontsize+2)
