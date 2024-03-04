
addpath ../colormaps/
load MAVS2_Phase_Average.mat

fontsize = 17;

figure(1)
set(gcf,'color','w')
subplot(2,1,1)
pcolor(phase,hab(:,1),temp)
shading flat;colorbar;colormap(redblue)
xlabel('Phase');ylabel('HAB (m)')
set(gca,'Fontsize',fontsize)
title('Temperature','Fontsize',fontsize+2)

subplot(2,1,2)
pcolor(phase,hab(:,1),N2)
shading flat;colorbar;colormap(redblue)
xlabel('Phase');ylabel('HAB (m)')
set(gca,'Fontsize',fontsize)
title('N^2','Fontsize',fontsize+2)
