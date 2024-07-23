%%% Calculate linear-fit shear, N2 and Ri of the large-scale flow

clear;

%--- MAVS1
load('MAVS1_shear.mat') % microseconds==>hours since 2021-07-06 04:45:32.999788
load('MAVS1_N2.mat')    % seconds since 2021-08-01 00:00:00
% load('MAVS1_shear_100m.mat') 
% load('MAVS1_N2_100m.mat')   
n2start = 1;
ustart = (31-6+1)*24*4-(4*4+3);

% %--- MAVS2
% % load('MAVS2_shear.mat') % microseconds==>hours since 2021-07-07 06:00:33.000323
% % load('MAVS2_N2.mat')    % seconds since 2021-07-06 00:00:00
% load('MAVS2_shear_100m.mat')
% load('MAVS2_N2_100m.mat')    
% n2start = 86400+6*3600+33;
% ustart = 1;

time_n2 = 1:length(time_temp); 
time_n2 = time_n2/3600; %%% convert to hours
n2idx = n2start:length(time_n2);
uend = ustart+time_n2(end)*4+1;
uidx = ustart:uend;

N2_zavg = N2_zavg(n2idx)';
smooth_N2_zavg = smooth_N2_zavg(n2idx)';
time_n2 = time_n2(n2idx);
time_n2 = time_n2 - time_n2(1);
shear_zavg = shear_zavg(uidx)';
time_uw = time_uw(uidx);
time_uw = time_uw-time_uw(1);

n2 = N2_zavg;
smooth_n2 = smooth_N2_zavg;
shear = shear_zavg;
time = time_n2;

%%% Interpolate 
shear_int = interp1(time_uw,shear,time_n2,'linear','extrap');

Ri = n2./(shear_int.^2);
smooth_Ri = smooth_n2./(shear_int.^2);


fontsize = 16;

figure(1)
clf;set(gcf,'Color','w')
subplot(3,1,1)
plot(time_n2/24,n2);
hold on;
plot(time_n2/24,smooth_n2,'--');
grid on;grid minor
axis tight
set(gca,'FontSize',fontsize);
title('N^2')
legend('Unsmoothed','Smoothed')
ylabel('(1/s^2)')
ylim([-2 13]*1e-6)

subplot(3,1,2)
% plot(time_uw/24,shear)
% hold on;
plot(time_n2/24,shear_int)
grid on;grid minor
axis tight
set(gca,'FontSize',fontsize);
title('Linear-fit shear')
ylabel('(1/s)')
ylim([-2.1 2.1]*1e-3)

subplot(3,1,3)
plot(time/24, 1./Ri);grid on;
hold on;
plot(time/24, 1./smooth_Ri,'--');
set(gca,'FontSize',fontsize);
legend('Unsmoothed','Smoothed')
xlabel('Time (days)')
grid on;grid minor
title('Inverse Richardson number')
axis tight
ylim([-2 4])


% figure(4)
% clf;
% plot(time,n2./max(abs(n2)));
% hold on;
% plot(time,shear_int./max(abs(shear_int)));
% plot(time, 1./Ri);
% % plot(time, 1./smooth_Ri,'--');
% ylim([-4 4])
% grid on;grid minor
% hold off;
% axis tight


clear ustart uend fontsize time_n2 time_uw shear time_temp n2start n2idx uidx N2_zavg shear_zavg smooth_N2_zavg

% save('MAVS1_Ri.mat')


