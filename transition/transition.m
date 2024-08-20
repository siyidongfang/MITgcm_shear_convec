clear;close all;
omega = 2*pi/43200;
N = 1e-3;
% shear_all = [0.1:0.1:2]*1e-3;
% Ri_all = N^2./shear_all.^2;
% Ri_all = 0.25:0.01:2;
inverseRi = 0:0.1:4;
Ri_all = 1./inverseRi;

Nt = 1;
dt = 5;
tt = dt:dt:43200*Nt;
Ln = length(Ri_all);
intA = NaN*zeros(1,Ln);
maxA = NaN*zeros(1,Ln);
k0=1;
rw = 0;
m0=rw*k0;

figure(1)
l1=plot(tt/3600,omega^2/2*ones(1,length(tt)),'k--','LineWidth',2);
hold on;
l2=plot(tt/3600,omega^2*ones(1,length(tt)),'b--','LineWidth',2);

% for n=1:21
for n=22:Ln
    Ri = Ri_all(n)
    B = abs(m0/k0*omega/N-Ri^(-0.5)*sin(omega*tt));
    C = 1./(omega^2/N^2+B.^2)-1;
    A = omega^2.*(1+C);

    intA(n) = sum(A.*dt);
    maxA(n) = max(A);
    minA(n) = min(A);
    meanA(n)=mean(A);

    intC(n) = sum(C.*dt);

    figure(1)
    plot(tt/3600,A);
    set(gca,'Fontsize',16)
    xlabel('Time (hours)')
    title('A(t), m0/k0=0')
    set(gca,'YScale', 'log')
    % ylim([0 1e-6])
    % ylim([-2 2])
end

legend([l1 l2],'\omega^2/2','\omega^2')


figure(3)
plot(inverseRi,minA)
hold on;
plot(inverseRi,omega^2/2*ones(1,length(inverseRi)))
set(gca,'YScale','log')
grid on;
grid minor;
set(gca,'fontsize',15)


aaa = meanA./omega^2;

figure(4)
plot(inverseRi,aaa)
set(gca,'YScale','log')
grid on;
grid minor;
