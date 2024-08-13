% close all;
omega = 2*pi/43200;
N = 1e-3;
shear_all = [0.1:0.1:2]*1e-3;
Ri_all = N^2./shear_all.^2;
% Ri_all = 0.25:0.01:2;

Nt = 1;
dt = 5;
tt = dt:dt:43200*Nt;
Ln = length(Ri_all);
intA = NaN*zeros(1,Ln);
maxA = NaN*zeros(1,Ln);
k0=1;
rw = 100;
m0=rw*k0;

for n=1:Ln
    Ri = Ri_all(n)
    A = omega^2./(omega^2/N^2+(m0/k0*omega/N-Ri^(-0.5)*sin(omega*tt)).^2);
    intA(n) = sum(A.*dt);
    maxA(n) = max(A);
    figure(4)
    hold on;
    plot(tt/3600,A);
    set(gca,'Fontsize',16)
    xlabel('Time (hours)')
    title('A(t), m0/k0=100')
    set(gca,'YScale', 'log')
end


% figure(2)
% plot(1./Ri_all,maxA);grid on;grid minor
% 
% figure(3)
% plot(1./Ri_all,intA);grid on;grid minor
% set(gca,'fontsize',15)
% ylim([0 0.1])
