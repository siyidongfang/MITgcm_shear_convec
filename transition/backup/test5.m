clear;close all;
Ptide = 43200;
omega=2*pi/43200;
Nt = 5;
dt = 5;
tt = dt:dt:43200*Nt;

y0 = sin(omega*tt).^2;

omega1 = 10*omega;
omega2 = 0.97*omega;


t1 =  dt:dt:0.125*Ptide;
t2 =  dt+0.125*Ptide:dt:0.375*Ptide;
t3 =  dt+0.375*Ptide:dt:0.625*Ptide;
t4 =  dt+0.625*Ptide:dt:0.875*Ptide;
t5 =  dt+0.875*Ptide:dt:Ptide;

y1 = sin(omega1*t1).^2;
init = y1(end)-sin(omega2*0.125*Ptide).^2;
y2 = init+sin(omega2*t2).^2;
init = y2(end)-sin(omega1*0.375*Ptide).^2;
y3 = init+sin(omega1*t3).^2;
init = y3(end)-sin(omega2*0.625*Ptide).^2;
y4 = init+sin(omega2*t4).^2;
init = y4(end)-sin(omega1*0.875*Ptide).^2;
y5 = init+sin(omega1*t5).^2;


y = [y1 y2 y3 y4 y5];

for n=1:Nt-1
    t1 = n*Ptide+ (dt:dt:0.125*Ptide);
    t2 = n*Ptide+ (dt+0.125*Ptide:dt:0.375*Ptide);
    t3 = n*Ptide+ (dt+0.375*Ptide:dt:0.625*Ptide);
    t4 = n*Ptide+ (dt+0.625*Ptide:dt:0.875*Ptide);
    t5 = n*Ptide+ (dt+0.875*Ptide:dt:Ptide);

    init = y(end)-sin(omega1*Ptide*n).^2;
    y1 = init+sin(omega1*t1).^2;
    init = y1(end)-sin(omega2*(0.125+n)*Ptide).^2;
    y2 = init+sin(omega2*t2).^2;
    init = y2(end)-sin(omega1*(0.375+n)*Ptide).^2;
    y3 = init+sin(omega1*t3).^2;
    init = y3(end)-sin(omega2*(0.625+n)*Ptide).^2;
    y4 = init+sin(omega2*t4).^2;
    init = y4(end)-sin(omega1*(0.875+n)*Ptide).^2;
    y5 = init+sin(omega1*t5).^2;

    y = [y y1 y2 y3 y4 y5];
end


figure(1)
plot(tt/3600/12,y0,'LineWidth',2)
hold on;
plot(tt/3600/12,y,'LineWidth',2)
set(gca,'Fontsize',18)
xlabel('Time (tidal cycles)')
grid on;grid minor;


