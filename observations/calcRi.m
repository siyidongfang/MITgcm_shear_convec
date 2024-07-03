%%% Calculate linear-fit shear, N2 and Ri of the large-scale flow

clear;

%--- MAVS1
load('MAVS1_shear.mat') % microseconds==>hours since 2021-07-06 04:45:32.999788
load('MAVS1_N2.mat')    % seconds since 2021-08-01 00:00:00
n2start = 1;
ustart = (31-6+1)*24*4-(4*4+3);

% %--- MAVS2
% load('MAVS2_shear.mat') % microseconds==>hours since 2021-07-07 06:00:33.000323
% load('MAVS2_N2.mat')    % seconds since 2021-07-06 00:00:00
% n2start = 86400+6*3600+33;
% ustart = 1;

time_n2 = 1:length(time_temp); 
time_n2 = time_n2/3600; %%% convert to hours
n2idx = n2start:length(time_n2);
uend = ustart+time_n2(end)*4+1;
uidx = ustart:uend;

N2_zavg = N2_zavg(n2idx)';
time_n2 = time_n2(n2idx);
time_n2 = time_n2 - time_n2(1);
shear_zavg = shear_zavg(uidx)';
time_uw = time_uw(uidx);
time_uw = time_uw-time_uw(1);

n2 = N2_zavg;
shear = shear_zavg;
time = time_n2;

%%% Interpolate 
shear_int = interp1(time_uw,shear,time_n2,'linear','extrap');

figure(1)
clf;plot(time_n2,n2)
grid on;grid minor
axis tight

figure(2)
clf;plot(time_uw,shear)
hold on;
plot(time_n2,shear_int)
grid on;grid minor
axis tight

clear time_n2 time_uw shear time_temp n2start n2idx uidx N2_zavg N2_zavg shear_zavg

Ri = n2./(shear_int.^2);

figure(3)
clf;
plot(time, 1./Ri);grid on;
ylim([0 1])

save('MAVS1_Ri.mat')

figure(4)
clf;
plot(time,n2./max(abs(n2)));
hold on;
plot(time,shear_int./max(abs(shear_int)));
plot(time, 1./Ri);
ylim([-10 10])
grid on;grid minor
hold off;


