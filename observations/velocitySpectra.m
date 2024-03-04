%%%
%%% velocitySpectra.m
%%%
%%% Plot temporal spectra of velocity and shear of the mooring data

clear;close all;

addpath /Users/ysi/MITgcm_shear_convec/observations/topography/
addpath /Users/ysi/MITgcm_shear_convec/analysis/colormaps/
addpath /Users/ysi/MITgcm_shear_convec/analysis/colormaps/cmocean/
addpath /Users/ysi/MITgcm_shear_convec/observations/moorings_gridded_adcp_v2/

%%% Calculate along-canyon velocity uu and cross-slope velocity vv
load("topog1D_28km.mat",'lat9','lon9','depth9')


ncfname = 'MAVS2_24606.nc';
lat = 54.18229903565569; 
lat1 = lat9(9);
lon1 = lon9(9);
lat2 = lat9(10);
lon2 = lon9(10);
t_idx = 32:8574;

% ncfname = 'MAVS1_24608.nc';
% lat = 54.197478060955085; 
% lat1 = lat9(8);
% lon1 = lon9(8);
% lat2 = lat9(9);
% lon2 = lon9(9);
% t_idx = 40:8756;


% ncfname = 'MP1_24839.nc';
% lat = 54.23853400073004;
% lat1 = lat9(5);
% lon1 = lon9(5);
% lat2 = lat9(6);
% lon2 = lon9(6);
% t_idx = 12:644;

tt = ncread(ncfname,'temperature')';
u_meridional = ncread(ncfname,'u')';
v_zonal = ncread(ncfname,'v')';
ww = ncread(ncfname,'w')'; 
pp = ncread(ncfname,'pressure')'; 
depth = ncread(ncfname,'depth')'; 
time  = double(ncread(ncfname,'time'))'/1e6; %%% convert microsecond to second

dy_canyon = distance(lat1,lon2,lat2,lon2,referenceEllipsoid('GRS80','km'));
dx_canyon = distance(lat2,lon1,lat2,lon2,referenceEllipsoid('GRS80','km'));

angle_canyon = atand(dy_canyon/dx_canyon); %%% The angle between the zonal direction and the canyon
uu = u_meridional*cosd(angle_canyon)-v_zonal*sind(angle_canyon);
vv = u_meridional*sind(angle_canyon)+v_zonal*cosd(angle_canyon);

%%% Calculate the kinetic energy
ke = 0.5*(u_meridional.^2+v_zonal.^2+ww.^2);

fontsize = 16;
t_idx_plot = 97:96*10;
% t_idx_plot = 1:length(time);

dz = diff(depth); %%% DOUBLE CHECK dz!!
dz = dz(1)
dudz = diff(uu)./dz;
% dvdz = diff(vv)./dz;
% dwdz = diff(ww)./dz;

save('MAVS2_velocity.mat','uu','vv','ww','u_meridional','v_zonal','depth','angle_canyon','time')



figure(1);clf;set(gcf,'color','w');
subplot(2,1,1)
plot(time/86400,pp);ylim([1150 1170]);grid on;grid minor;
title('Pressure (dbar)');xlabel('time (days)')
set(gca,'FontSize',fontsize)
subplot(2,1,2)
plot(time/86400,tt);grid on;grid minor;
title('Temperature (^oC)');xlabel('time (days)')
set(gca,'FontSize',fontsize)

figure(2);clf;set(gcf,'color','w');
subplot(3,1,1)
pcolor(time(t_idx_plot)/86400,depth,u_meridional(:,t_idx_plot));
shading interp;colorbar;colormap(redblue);clim([-0.4 0.4])
title('Zonal velocity (m/s)')
ylabel('Depth (m)');xlabel('time (days)');set(gca,'FontSize',fontsize)
subplot(3,1,2)
pcolor(time(t_idx_plot)/86400,depth,v_zonal(:,t_idx_plot));
shading interp;colorbar;colormap(redblue);clim([-0.4 0.4])
title('Meridional velocity (m/s)')
ylabel('Depth (m)');xlabel('time (days)');set(gca,'FontSize',fontsize)
subplot(3,1,3)
pcolor(time(t_idx_plot)/86400,depth,ww(:,t_idx_plot));
shading interp;colorbar;colormap(redblue);clim([-0.1 0.1])
title('Vertical velocity (m/s)')
ylabel('Depth (m)');xlabel('time (days)');set(gca,'FontSize',fontsize)


figure(3);clf;set(gcf,'color','w');
subplot(3,1,1)
pcolor(time(t_idx_plot)/86400,depth,uu(:,t_idx_plot));
shading interp;colorbar;colormap(redblue);clim([-0.4 0.4])
title('Up-canyon velocity (m/s)')
ylabel('Depth (m)');xlabel('time (days)');set(gca,'FontSize',fontsize)

subplot(3,1,2)
pcolor(time(t_idx_plot)/86400,depth,vv(:,t_idx_plot));
shading interp;colorbar;colormap(redblue);clim([-0.4 0.4])
title('Cross-canyon velocity (m/s)')
ylabel('Depth (m)');xlabel('time (days)');set(gca,'FontSize',fontsize)

subplot(3,1,3)
pcolor(time(t_idx_plot)/86400,depth,ww(:,t_idx_plot));
shading interp;colorbar;colormap(redblue);clim([-0.1 0.1])
title('Vertical velocity (m/s)')
ylabel('Depth (m)');xlabel('time (days)');set(gca,'FontSize',fontsize)

figure(4)
pcolor(time(t_idx_plot)/86400,depth,ke(:,t_idx_plot));
shading interp;colorbar;colormap(redblue);clim([-0.1 0.1])
title('Kinetic Energy (m^2/s^2)')
ylabel('Depth (m)');xlabel('time (days)');set(gca,'FontSize',fontsize)


%%

%%% Length of record
lev = 15;
time_ftt = time(t_idx)/86400;
Ntime = length(time_ftt);

%%% Vectors of spectral period and frequency
period = diff(time_ftt([1 end])) ./ (0:1:Ntime/2-1);
freq = 2*pi./period;

%%% Coriolis parameter
f = 2*2*pi/86164*sind(lat) *86400;
period_m2 = 12.4206012*3600;
period_k1 = 23.93447213*3600;
freq_m2 = 2*pi/period_m2 *86400;
freq_k1 = 2*pi/period_k1 *86400;
 
aaa = uu(lev,t_idx);
sum(isnan(aaa))

%%%%%%%%%
%%%%%%%%%
%%%%%%%%% TO DO: EXCLUDE NANs FROM VELOCITY/SHEAR
%%%%%%%%%
%%%%%%%%%

%%% Compute FFTs
tfft = fft(tt(t_idx));
% sfft = fft(A253.mat.s);
ufft = fft(uu(lev,t_idx));
vfft = fft(vv(lev,t_idx));
ushearfft = fft(dudz(lev,t_idx));
freq_ref =  2;

%%% Make plots
% figure(5)
% clf;set(gcf,'color','w');
% loglog(period,abs(tfft(1:Ntime/2)).^2);
% hold on;
% loglog(abs(2*pi./[f f]),[1e-3 1e7],'k:','LineWidth',2);
% loglog(abs(2*pi./[freq_m2 freq_m2]),[1e-3 1e7],'r:','LineWidth',2);
% loglog(abs(2*pi./[freq_k1 freq_k1]),[1e-3 1e7],'b:','LineWidth',2);
% hold off;
% xlabel('Period (days)');
% ylabel('Spectral power (^oC^2)');
% grid on;grid minor;set(gca,'FontSize',fontsize)
% ylim([1e-3 2e6])
% title('Temporal spectra of temperature')


figure(6)
clf;set(gcf,'color','w');
loglog(freq,abs(tfft(1:Ntime/2)).^2);
hold on
loglog(abs([f f]),[1e-3 1e7],'k:','LineWidth',2);
loglog(abs([freq_m2 freq_m2]),[1e-3 1e7],'r:','LineWidth',2);
loglog(abs([freq_k1 freq_k1]),[1e-3 1e7],'b:','LineWidth',2);
loglog(freq,5e4*(freq/freq_ref).^(-2),'k--');
hold off;
xlabel('Frequency (rad/day)');
ylabel('Spectral power (^oC^2)');
grid on;grid minor;set(gca,'FontSize',fontsize)
% ylim([1e-3 2e6])
title('Temporal spectra of temperature')

%%%!!!!!!!!!!DOUBLE-CHECK KE SPECTRA CALCULATION
% figure(7)
% clf;set(gcf,'color','w');
% loglog(period,abs(ufft(1:Ntime/2)).^2+abs(vfft(1:Ntime/2)).^2)
% % loglog(period,abs(ufft(1:Ntime/2)).^2)
% hold on
% loglog(abs(2*pi./[f f]),[1e-4 1e7],'k:','LineWidth',2);
% loglog(abs(2*pi./[freq_m2 freq_m2]),[1e-3 1e7],'r:','LineWidth',2);
% loglog(abs(2*pi./[freq_k1 freq_k1]),[1e-3 1e7],'b:','LineWidth',2);
% hold off;
% xlabel('Period (days)');
% ylabel('Spectral power (m^2/s^2)');
% grid on;grid minor;set(gca,'FontSize',fontsize)
% ylim([1e-1 5e5])
% title('Temporal spectra of kinetic energy')


% %%%!!!!!!!!!!DOUBLE-CHECK KE SPECTRA CALCULATION
% figure(8)
% clf;set(gcf,'color','w');
% loglog(freq,abs(ufft(1:Ntime/2)).^2+abs(vfft(1:Ntime/2)).^2)
% hold on
% loglog(abs([f f]),[1e-3 1e7],'k:','LineWidth',2);
% loglog(abs([freq_m2 freq_m2]),[1e-3 1e7],'r:','LineWidth',2);
% loglog(abs([freq_k1 freq_k1]),[1e-3 1e7],'b:','LineWidth',2);
% loglog(freq,3e4*(freq/freq_ref).^(-2),'k--')
% hold off;
% xlabel('Frequency (rad/day)');
% ylabel('Spectral power (m^2/s^2)');
% grid on;grid minor;set(gca,'FontSize',fontsize)
% ylim([1e-1 5e5])
% title('Temporal spectra of kinetic energy????')


figure(9)
clf;set(gcf,'color','w');
loglog(freq,abs(ufft(1:Ntime/2)).^2)
hold on
loglog(freq,abs(vfft(1:Ntime/2)).^2)
loglog(abs([f f]),[1e-3 1],'k:','LineWidth',2);
loglog(abs([freq_m2 freq_m2]),[1e-3 1],'r:','LineWidth',2);
loglog(abs([freq_k1 freq_k1]),[1e-3 1],'b:','LineWidth',2);
% loglog(freq,3e4*(freq/freq_ref).^(-2),'k--')
hold off;
xlabel('Frequency (rad/day)');
ylabel('Spectral power (m^2/s^2)');
grid on;grid minor;set(gca,'FontSize',fontsize)
% ylim([0e-1 5e5])
ylim([0e-1 5e2])
title('Temporal spectra of kinetic energy')
legend('Along-canyon','Cross-canyon')


figure(10)
clf;set(gcf,'color','w');
loglog(freq,abs(ushearfft(1:Ntime/2)).^2)
hold on
loglog(abs([f f]),[1e-6 0.1],'k:','LineWidth',2);
loglog(abs([freq_m2 freq_m2]),[1e-6 0.1],'r:','LineWidth',2);
loglog(abs([freq_k1 freq_k1]),[1e-6 0.1],'b:','LineWidth',2);
% loglog(freq,3e4*(freq/freq_ref).^(-2),'k--')
hold off;
xlabel('Frequency (rad/day)');
ylabel('Spectral power (1/s^2)');
grid on;grid minor;set(gca,'FontSize',fontsize)
ylim([1e-6 10])
title('Temporal spectra of along-canyon velocity shear')
xlim([0.05 500])




