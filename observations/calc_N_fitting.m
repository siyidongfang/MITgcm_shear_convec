clear all;close all

% fname = '_bottom_begin'
fname = '_full_begin'

%--- MAVS2
load('MAVS2_LinearShear.mat') %% Linear-fit shear
load('MAVS2_N2.mat')    % seconds since 2021-07-06 00:00:00
% load('MAVS2_shear.mat') % microseconds==>hours since 2021-07-07 06:00:33.000323

% load('MAVS2_LinearShear_100m.mat') %% Linear-fit shear
% load('MAVS2_N2_100m.mat')  
% % load('MAVS2_shear_100m.mat')

shear_zavg = shear_linear;
time_uw=time/3600;
n2start = 86400+6*3600+33;
ustart = 1;
fontsize = 19;

load_colors;

time_n2 = 1:length(time_temp); 
time_n2 = time_n2/3600; %%% convert to hours
n2idx = n2start:length(time_n2);
uend = ustart+time_n2(end)*4+1;
uidx = ustart:uend;

N2_zavg = N2_zavg(n2idx)';
smooth_N2_zavg = smooth_N2_zavg(n2idx)';
time_n2 = time_n2(n2idx);
time_n2 = time_n2 - time_n2(1);
shear_zavg = shear_zavg(uidx)';
time_uw = time_uw(uidx);
time_uw = time_uw-time_uw(1);

n2 = N2_zavg;
smooth_n2 = smooth_N2_zavg;
shear = shear_zavg;
time = time_n2;


%%% Interpolate 
shear_int = interp1(time_uw,shear,time_n2,'linear','extrap');

%%% Curve fitting to a sinusoidal function
Period_m2 = 12*3600+25*60;
omega_m2 = 2*pi/Period_m2;
% phi = 0;
% Lambda_fit = 1e-3;
% shear_fit = Lambda_fit*cos(omega_m2*time_n2*3600+phi);

x = time_n2'*3600/86400;
y = shear_int';
% Nstarty = find(~isnan(y),1);
% Nstarty=round(length(y)/4*3)-3600-86400;
% Nendy=length(y)-43200;

Nstarty = 86400/2;
Nendy=round(length(y)/4);
yidx = Nstarty:Nendy;
y=y(yidx);
x=x(yidx);
fo = fitoptions('Method','NonlinearLeastSquares');
mdl = fittype('a*sin(b*x+c)+d','indep','x','options',fo);
fittedmdl_shear= fit(x,y,mdl,'start',[rand(),omega_m2*86400,-pi/2,rand()])
% mdl = fittype('a*sin(b*x+c)','indep','x','options',fo);
% fittedmdl_shear = fit(x,y,mdl,'start',[rand(),omega_m2*86400,-pi/2])

xfit_shear = x;
yfit_shear = fittedmdl_shear(xfit_shear);


figure(1);
clf;set(gcf,'Color','w')
hold on;
plot(x,y,'Color','k','LineWidth',1.5)
% l1=plot(fittedmdl_shear);
plot(xfit_shear,yfit_shear,'Color',yellow,'LineWidth',2);
grid on;grid minor
axis tight
set(gca,'FontSize',fontsize);
title('Linear-fit shear')
ylabel('(1/s)')
ylim([-2.5 2.5]*1e-3)
xlabel('Time (days)')
legend('Linear-fit shear','Sinusoidal curve fitting')
box on;
% xlim([10 14.2])
xlim([0.5 3.7])
% print('-dpng','-r200',['figures/fig1' fname '.png']);



y = smooth_n2';
y=y(yidx);
mdl = fittype('a*sin(b*x+c)+d','indep','x','options',fo);
fittedmdl_N2 = fit(x,y,mdl,'start',[rand(),omega_m2*86400,-pi/4,rand()])
xfit_N2 = x;
yfit_N2 = fittedmdl_N2(xfit_N2);


Ri_curvefit = yfit_N2./(yfit_shear.^2);


figure(2);
clf;set(gcf,'Color','w')
hold on;
plot(time_n2(yidx)/24,n2(yidx),'Color',gray,'LineWidth',1.5)
plot(time_n2(yidx)/24,smooth_n2(yidx),'--','Color','k','LineWidth',1.5)
plot(xfit_N2,yfit_N2,'Color',yellow,'LineWidth',2);
grid on;grid minor
axis tight
set(gca,'FontSize',fontsize);
title('N^2')
legend('Unsmoothed N^2','Smoothed N^2','Sinusoidal curve fitting')
ylabel('(1/s^2)')
ylim([-1 13]*1e-6)
xlabel('Time (days)')
box on;
% xlim([10 14.2])
xlim([0.5 3.7])
% print('-dpng','-r200',['figures/fig2' fname '.png']);






Ri = n2./(shear_int.^2);
smooth_Ri = smooth_n2./(shear_int.^2);

%%%%%% Day 0.5-4.2
% %%% full depth
% shear_fit = 0.0008605;
% N2 = 4.316e-06;
% % N2 = 2.917e-06;
%%% Bottom 100m
shear_fit = 0.0005052;
N2 = 3.263e-06 ;
% N2 = 2.237e-06;

% %%%%%% Day 10-14.2
% % %%% full depth 
% % shear_fit = 0.0004248;
% % N2 = 3.422e-06;
% % % N2 = 1.335e-06;
% 
% % %%% Bottom 100m
% % shear_fit = 0.0008715;
% % N2 = 3.571e-06 ;
% % % N2 = 2.467e-06;

dt = 600;
NTtide = 100;
omega = omega_m2;
Nt = NTtide/omega/dt;

shear_convec = cosd(topo)/sind(topo)*omega;
dt_ri = dt/1000;
tt_ri = dt_ri:dt_ri:Nt*dt;

shear = shear_fit;
Ri_inverse = (shear*cos(omega*tt_ri)).^2./(N2*cosd(topo) - N2*sind(topo)/omega*shear*sin(omega*tt_ri));
Ri_min = 1/max(Ri_inverse) 


figure(3)
clf;set(gcf,'Color','w')
plot(time(yidx)/24, 1./Ri(yidx),'Color',gray,'LineWidth',1.5);grid on;
hold on;
plot(time(yidx)/24, 1./smooth_Ri(yidx),'--','Color','k','LineWidth',1.5);
plot(time(yidx)/24, 1./Ri_min*ones(1,length(yidx)),'-','LineWidth',2);
plot(time(yidx)/24, 1./Ri_curvefit,'--','LineWidth',2);
set(gca,'FontSize',fontsize);
legend('Unsmoothed','Smoothed')
xlabel('Time (days)')
grid on;grid minor
title('Inverse Richardson number')
axis tight
ylim([-0.5 5])
% xlim([10 14.2])
xlim([0.5 3.7])
legend('Inverse Ri (unsmoothed N^2)', 'Inverse Ri (smoothed N^2)','Inverse minimum Ri, sinusoidal curve fitting')

% print('-dpng','-r200',['figures/fig3' fname '.png']);


clear ustart uend fontsize time_n2 time_uw shear time_temp n2start n2idx uidx N2_zavg shear_zavg smooth_N2_zavg



