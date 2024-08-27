
clear;close all;

omega = 2*pi/43200;
N = 1e-3;

inverseRi = 0:0.1:4;

A = N^2./(1+N^2/omega^2*inverseRi);

plot(inverseRi,A,'LineWidth',2)
hold on;
l2 =plot(inverseRi,omega^2/2*ones(1,length(inverseRi)),'LineWidth',2);
set(gca,'YScale','log')
grid on;
grid minor;
set(gca,'Fontsize',17)
title('Minimum A(t) as a function of R_i^{-1} (m_0/k_0~0)')
xlabel('R_i^{-1}')
ylabel('(s^{-2})')
legend(l2,'\omega^2/2')
