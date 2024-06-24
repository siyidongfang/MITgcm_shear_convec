
addpath ../analysis/colormaps/
addpath ../observations/topography/

load_colors;

ncfname = 'blt_canyon_mb_qc.nc';
lat = ncread(ncfname,'lat')';
lon = ncread(ncfname,'lon')';
z = ncread(ncfname,'z');

fontsize = 16;

figure(15);
clf;set(gcf,'color','w');
hold on;
surf(lon,lat,-z','EdgeColor','None');
shading flat; 
%colorbar;
colormap(flip(cmocean('rain')))
xlabel('Longitude');
ylabel('Latitude');
% zlabel('z (km)','Rotation',0);
title('Canyon');
set(gca,'FontSize',fontsize)
box on;grid on;grid minor;
view([-77.9824 57.7483]);
axis tight;
clim([-2500 -500])
zlim([-2500 -500])
% print('-djpeg','-r200',[figdir 'canyon.jpeg']);


% pbaspect([Lx/Ly 1 1]);
