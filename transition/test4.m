clear;close all;
Ptide = 43200;
omega=2*pi/43200;
Nt = 5;
dt = 5;
tt = dt:dt:43200*Nt;

y1 = sin(omega*tt).^2;

omega1 = 3*omega;
omega2 = 0.9*omega;

alpha = 0.6;

t21 =  dt:dt:alpha*Ptide;
t22 =  dt+alpha*Ptide:dt:Ptide;

y21 = sin(omega1*t21).^2;
init = y21(end)-sin(omega2*alpha*Ptide).^2;
y22 = init+sin(omega2*t22).^2;
y2 = [y21 y22];

for n=1:Nt-1
    t23 =  n*Ptide+t21;
    t24 =  n*Ptide+t22;
    init = y2(end)-sin(omega1*Ptide*n).^2;
    y23 = init+sin(omega1*t23).^2;
    init = y23(end)-sin(omega2*(alpha+n)*Ptide).^2;
    y24 = init+sin(omega2*t24).^2;
    y2 = [y2 y23 y24];
end


figure(1)
plot(tt/3600/12,y1,'LineWidth',2)
hold on;
plot(tt/3600/12,y2,'LineWidth',2)
set(gca,'Fontsize',18)
xlabel('Time (tidal cycles)')
grid on;grid minor;


