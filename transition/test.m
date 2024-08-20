clear;
close all;
a=2;
b=3;

Nt = 1;
dt = 5;
tt = dt:dt:43200*Nt;
omega = 2*pi/43200;

% B = abs(a-b*sin(omega*tt));
B = 1./(1-cos(2*omega*tt));

figure(4)
plot(tt/3600,B);
set(gca,'Fontsize',16)
xlabel('Time (hours)')
% ylim([0 1]*1e-6)