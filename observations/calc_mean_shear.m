%%%
%%% calc_mean_shear.m
%%%
%%% Calculate the mean shear from observations

load('MAVS2_velocity.mat')

nStart = 263;
nEnd = 454;
tidx = nStart:nEnd;
zidx = 4:18;
uselect = uu(zidx,tidx)';
vselect = vv(zidx,tidx)';
wselect = ww(zidx,tidx)';
time_uw = time(tidx)/3600; %%% in hours
time_uw = time_uw-time_uw(1);
depth_uw = depth(zidx);

fontsize = 20;

shear = diff(uselect,1,2)./diff(-depth_uw);
depth_shear = 0.5*(depth_uw(1:end-1)+depth_uw(2:end));


figure(4);
clf;set(gcf,'Color','w');set(gcf,'Position', [56 352 600 305]);
pcolor(time_uw,depth_uw,uselect');shading flat;colorbar;colormap(redblue)
xlabel('Time (hours)');ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);axis ij;clim([-0.5 0.5])
title('u (m/s)')

figure(5);
clf;set(gcf,'Color','w');set(gcf,'Position', [56 352 600 305]);
pcolor(time_uw,depth_shear,shear');shading flat;colorbar;colormap(redblue)
xlabel('Time (hours)');ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);axis ij;clim([-5 5]/1000)
title('Shear (1/s)')

%%
%%% Apply a low-pass filter to velocity, and then calculate the mean shear
% u_smooth = zeros(size(uselect));
% span = 10;
% for nt = 1:length(uselect)
%     u_smooth(nt,:) = smooth(uselect(nt,:),span)';
% end

max_shear_nosmooth = max(abs(shear));

windowsize = round(100./mean(diff(depth_uw)));
u_smooth = smoothdata(uselect,2,'gaussian',windowsize,'omitnan');
shear_smooth = diff(u_smooth,1,2)./diff(-depth_uw);
max_shear_100 = max(abs(shear_smooth));

windowsize = round(200./mean(diff(depth_uw)));
u_smooth = smoothdata(uselect,2,'gaussian',windowsize,'omitnan');
shear_smooth = diff(u_smooth,1,2)./diff(-depth_uw);
max_shear_200 = max(abs(shear_smooth));

windowsize = round(300./mean(diff(depth_uw)));
u_smooth = smoothdata(uselect,2,'gaussian',windowsize,'omitnan');
shear_smooth = diff(u_smooth,1,2)./diff(-depth_uw);
max_shear_300 = max(abs(shear_smooth));

windowsize = round(400./mean(diff(depth_uw)));
u_smooth = smoothdata(uselect,2,'gaussian',windowsize,'omitnan');
shear_smooth = diff(u_smooth,1,2)./diff(-depth_uw);
max_shear_400 = max(abs(shear_smooth));

figure(6);
clf;set(gcf,'Color','w');set(gcf,'Position', [56 352 600 305]);
pcolor(time_uw,depth_uw,uselect');shading flat;colorbar;colormap(redblue)
xlabel('Time (hours)');ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);axis ij;clim([-0.5 0.5])
title('Smoothed u (m/s)')

figure(7);
clf;set(gcf,'Color','w');set(gcf,'Position', [56 352 600 305]);
pcolor(time_uw,depth_shear,shear_smooth');shading flat;colorbar;colormap(redblue)
xlabel('Time (hours)');ylabel('Depth (m)')
set(gca,'Fontsize',fontsize);axis ij;clim([-5 5]/1000)
title('Smoothed shear (1/s)')


figure(8)
clf;set(gcf,'Color','w');
l0 = plot(max_shear_nosmooth,depth_shear,'LineWidth',2,'Color','k');
hold on;
l1=plot(max_shear_100,depth_shear,'LineWidth',2);
l2=plot(max_shear_200,depth_shear,'LineWidth',2);
l3=plot(max_shear_300,depth_shear,'LineWidth',2);
l4=plot(max_shear_400,depth_shear,'LineWidth',2);
grid on;grid minor;
set(gca,'Fontsize',fontsize);axis ij;
xlabel('Velocity shear (1/s)');ylabel('Depth (m)')
legend([l0 l1 l2 l3 l4],'Velocity shear (no smooth)','Window size = 100 m',...
    'Window size = 200 m','Window size = 300 m','Window size = 400 m')






