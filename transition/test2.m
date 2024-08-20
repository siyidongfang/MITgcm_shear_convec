
clear;close all;

omega = 2*pi/43200;
N = 1e-3;

inverseRi = 0:0.1:4;

A = N^2./(1+N^2/omega^2*inverseRi);

plot(inverseRi,A)
hold on;
plot(inverseRi,omega^2/2*ones(1,length(inverseRi)))
set(gca,'YScale','log')
grid on;
grid minor;
