clear;close all;
Ptide = 43200;
omega=2*pi/43200;
Nt = 10;
dt = 5;
tt = dt:dt:43200*Nt;

y1 = sin(omega*tt);

omega1 = 0.99*omega;
omega2 = 5*omega;

alpha = 0.2;


t21 =  dt:dt:alpha*Ptide;
t22 =  dt+alpha*Ptide:dt:Ptide;

y21 = sin(omega1*t21);
init = y21(end)-sin(omega2*alpha*Ptide);
y22 = init+sin(omega2*t22);
y2 = [y21 y22];

for n=1:Nt-1
    t23 =  dt+n*Ptide:dt:(alpha+n)*Ptide;
    t24 =  dt+(alpha+n)*Ptide:dt:(n+1)*Ptide;
    init = y2(end)-sin(omega1*Ptide*n);
    y23 = init+sin(omega1*t23);
    init = y23(end)-sin(omega2*(alpha+n)*Ptide);
    y24 = init+sin(omega2*t24);
    y2 = [y2 y23 y24];
end



figure(1)
plot(tt/3600/12,y1,'LineWidth',2)
hold on;
plot(tt/3600/12,y2,'LineWidth',2)
set(gca,'Fontsize',18)
xlabel('Time (tidal cycles)')
grid on;grid minor;


