%%% Calculate linear-fit shear, N2 and Ri of the large-scale flow

% %--- MAVS1
% load('MAVS1_shear.mat')
% load('MAVS1_N2.mat')

%--- MAVS2
load('MAVS2_shear.mat')
load('MAVS2_N2.mat')

time_n2 = 1:length(time_temp);
time_n2 = time_n2/86400; %%% convert to days

figure(1)
plot(time_n2,N2_zavg)
grid on;grid minor
axis tight

figure(2)
plot(time_uw,shear_zavg)
grid on;grid minor
axis tight

