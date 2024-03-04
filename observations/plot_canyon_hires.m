
clear;

addpath ../colormaps/
addpath ../observations/


load_colors;

ncfname = 'blt_canyon_mb_qc.nc';
lat = ncread(ncfname,'lat')';
lon = ncread(ncfname,'lon')';
z = ncread(ncfname,'z');

dy = diff(lat);dy = dy(1);
dx = diff(lon);dx = dx(1);

load Rockall_gebco.mat
lat_max = max(lat)+130*dy;
lat_min = min(lat)-170*dy;
lon_max = max(lon)+650*dx;
lon_min = min(lon)-130*dx; 

yy = [lat_min:dy:(lat(1)-dy) lat (lat(end)+dy):dy:lat_max];
xx = [lon_min:dx:(lon(1)-dx) lon (lon(end)+dx):dx:lon_max];

[YY,XX]=meshgrid(yy,xx);


lat_idx = find((ele_lat<=lat_max).*(ele_lat>=lat_min)==1);
lon_idx = find((ele_lon<=lon_max).*(ele_lon>=lon_min)==1);
lat_canyon = ele_lat(lat_idx);
lon_canyon = ele_lon(lon_idx);
elev_canyon = elev(lon_idx,lat_idx);

fontsize = 16;

figure(9)
contour(lon_canyon,lat_canyon,elev_canyon',[-3200:25:-200])
shading flat;colorbar;
shading flat;colorbar;
colormap(jet);
set(gca,'FontSize',fontsize);
title('Canyon topography (m): low resolution');
xlabel('Longitude'); ylabel('Latitude')
grid on;grid minor;



%%

load CTD_stations.mat
for nf = 1:15
    lat_CTD(nf) = CTD.Stn(nf).lat;
    lon_CTD(nf) = CTD.Stn(nf).lon;
end

%%% Make figures

%%% Extract topography of the canyon center based on the location of the 9 CTD casts.
lat9 = flip([lat_CTD(8:9) mean(lat_CTD(1:7)) lat_CTD(10:15)]);
lon9 = flip([lon_CTD(8:9) mean(lon_CTD(1:7)) lon_CTD(10:15)]);

% lat9 = [54.365 lat9 54.13 54.1];
% lon9 = [-12.075 lon9 -11.765 -11.765];
% lat9 = [54.396 54.378 54.365 lat9(1:2) 54.2764 lat9(3:end) 54.13];
% lon9 = [-12.13 -12.0898 -12.075 lon9(1:2) -11.9965 lon9(3:end) -11.765];
% lat9(9) = 54.226;
% lat9(11) = 54.203;

lat9 = [lat9(1:2) 54.2764 lat9(3:end) 54.13];
lon9 = [lon9(1:2) -11.9965 lon9(3:end) -11.765];
lat9(6) = 54.226;
lat9(8) = 54.203;

for i=1:length(lat9)-1
    lat1 = lat9(i);
    lon1 = lon9(i);
    lat2 = lat9(i+1);
    lon2 = lon9(i+1);
    [d1km(i),az(i)] = distance(lat1,lon1,lat2,lon2,referenceEllipsoid('GRS80','km'));
end

along_canyon = [0 cumsum(d1km)];

for i=1:length(lat9)
    lat_i = lat9(i);
    lon_i = lon9(i);
    [a,lat9idx(i)] = min(abs(lat-lat_i));
    [b,lon9idx(i)] = min(abs(lon-lon_i));
    depth9(i) = z(lon9idx(i),lat9idx(i));
end

%%% Linear interpolation for latitude and longitude along the canyon 
%%% Find the local minimum of depth as the center of the canyon
ll = 1:length(lat9);
llf = 1:1/100:length(lat9);
latn = interp1(ll,lat9,llf,'linear','extrap');
lonn = interp1(ll,lon9,llf,'linear','extrap');

ll = 1:length(lat);
llf = 1:1/10:length(lat);
latf = interp1(ll,lat,llf,'linear','extrap');
ll = 1:length(lon);
llf = 1:1/10:length(lon);
lonf = interp1(ll,lon,llf,'linear','extrap');
[Y,X] = meshgrid(lat,lon);
[Yq,Xq] = meshgrid(latf,lonf);
zf = interp2(Y,X,z,Yq,Xq,'linear');


for i=1:length(latn)
    lat_i = latn(i);
    lon_i = lonn(i);
    [a,latidx(i)] = min(abs(latf-lat_i));
    [b,lonidx(i)] = min(abs(lonf-lon_i));

    nx=0;
    % lonidx_max = lonidx(i)-nx:lonidx(i)+nx;
    % depth_max = z(lonidx_max,latidx(i));
    latidx_max = latidx(i)-nx:latidx(i)+nx;
    depth_max = zf(lonidx(i),latidx_max);
    [a,b] = max(depth_max);
    % lonidx(i) = lonidx_max(b);
    latidx(i) = latidx_max(b);
    depthn(i) = zf(lonidx(i),latidx(i));
end

latn = latf(latidx);
lonn = lonf(lonidx);

for i=1:length(latn)-1
    lat1 = latn(i);
    lon1 = lonn(i);
    lat2 = latn(i+1);
    lon2 = lonn(i+1);
    [d1kmn(i),az(i)] = distance(lat1,lon1,lat2,lon2,referenceEllipsoid('GRS80','km'));
end

along_canyonn = [0 cumsum(d1kmn)];



%%
figure(1);
clf;set(gcf,'color','w','Position',[237 206 592 360]);
set(gca,'Position',[0.12 0.14 0.7 0.75])
hold on;
% pcolor(lon,lat,z');
contour(lon,lat,z',[400:25:3000],'LineWidth',0.5);
contour(lon,lat,z',[400:100:3000],'LineWidth',2,'ShowText','on');
hold off;
shading flat;
handle = colorbar;
set(handle,'Position',[0.85 0.23 0.02 0.6], 'YDir', 'reverse');
colormap(cmocean('rain'));
set(gca,'FontSize',fontsize);
title('Canyon bathymetry (m)','FontSize',fontsize+5);
xlabel('Longitude'); ylabel('Latitude')
xlim([-12.13 -11.8])
hold on;scatter(lonn,latn,150,".",'LineWidth',1);
scatter(lon_CTD(8:15),lat_CTD(8:15),150,"x",'LineWidth',4);
scatter(lon_CTD(1:7),lat_CTD(1:7),150,"o",'LineWidth',4);
box on;grid on;grid minor;
% clim([500 2500])
% figdir = '/Users/csi/MITgcm_BLT/analysis/NCAR_proposal/';
% print('-djpeg','-r200',[figdir 'Bathymetry_canyon.jpeg']);


depthn = smooth(smooth(smooth(depthn)')')';
depth9_smooth = smooth(depth9)';
s_topog9_smooth = -diff(depth9_smooth)./d1km/1000;

s_topog9 = -diff(depth9)./d1km/1000;
along_mid = 0.5*(along_canyon(2:end)+along_canyon(1:end-1));

s_topogn = -diff(depthn)./d1kmn/1000;
along_midn = 0.5*(along_canyonn(2:end)+along_canyonn(1:end-1));

figure(2)
clf;set(gcf,'color','w');
subplot(2,1,1)
% plot(along_canyon,depth9);
% xlabel('Along-canyon distance (km)')
title('Canyon depth (m)')
hold on;
plot(along_canyonn,depthn);
axis ij;grid on;grid minor;
set(gca,'FontSize',fontsize);
subplot(2,1,2)
% plot(along_mid,s_topog9);
hold on;
plot(along_midn,abs(s_topogn));
xlabel('Along-canyon distance (km)')
title('Topographic slope')
grid on;grid minor;
ylim([0 0.33])
set(gca,'FontSize',fontsize);



mean_slope = (depth9(1)-depth9(9))/along_canyon(9)/1000;

mean_slope1 = (depth9(3)-depth9(5))/(along_canyon(5)-along_canyon(3))/1000;
mean_slope2 = (depth9(5)-depth9(9))/(along_canyon(9)-along_canyon(5))/1000;

Ttide = 44712;
omega_tides_obs = 2*pi/Ttide;
N_mean_obs = 0.62*2.4e-3;
fobs = 4*pi/86164*sind(54.2);
r_iw_obs = sqrt((omega_tides_obs^2-fobs^2)/(N_mean_obs^2-omega_tides_obs^2))


save('topography/topog1D_28km.mat','depth9','along_canyon','d1km','along_mid','s_topog9', ...
    'depthn','along_canyonn','d1kmn','along_midn','s_topogn',...
    'r_iw_obs','fobs','N_mean_obs','omega_tides_obs','mean_slope1','mean_slope2',...
    'lat9','lon9','depth9')


