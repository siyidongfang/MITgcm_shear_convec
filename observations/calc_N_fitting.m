clear all;close all

% fname = '_bottom_begin'
fname = '_full_begin'

%--- MAVS2
load('MAVS2_LinearShear.mat') %% Linear-fit shear % seconds==>hours since 2021-07-07 06:00:33.000323
load('MAVS2_N2.mat')    % seconds since 2021-07-06 00:00:00
% load('MAVS2_shear.mat') % microseconds==>hours since 2021-07-07 06:00:33.000323

% load('MAVS2_LinearShear_100m.mat') %% Linear-fit shear
% load('MAVS2_N2_100m.mat')  
% % load('MAVS2_shear_100m.mat')


nanidx = find(isnan(shear_linear(32:8574)));
shear_linear(nanidx+31) = 0.5* (shear_linear(nanidx+30) +  shear_linear(nanidx+32) );
% nanidx = find(isnan(shear_linear(32:8574)));
% shear_linear(nanidx+31) = 0.5* (shear_linear(nanidx+30) +  shear_linear(nanidx+32) );

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
shear_linear = shear_linear(uidx)';
time_uw = time_uw(uidx);
time_uw = time_uw-time_uw(1);

n2 = N2_zavg;
smooth_n2 = smooth_N2_zavg;
shear = shear_linear;
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

Nstarty = 86400/2+0.04*86400;
Nendy=round(length(y)/4);
yidx = Nstarty:Nendy;
y=y(yidx);
x=x(yidx);
x=x-x(1);
fo = fitoptions('Method','NonlinearLeastSquares');
mdl = fittype('a*sin(b*x+c)+d','indep','x','options',fo);
fittedmdl_shear= fit(x,y,mdl,'start',[rand(),omega_m2*86400,0,rand()])
% mdl = fittype('a*sin(b*x+c)','indep','x','options',fo);
% fittedmdl_shear = fit(x,y,mdl,'start',[rand(),omega_m2*86400,-pi/2])

xfit_shear = x;
yfit_shear = fittedmdl_shear(xfit_shear);

shear_fit = fittedmdl_shear.a;

shear_theory = shear_fit*sin(omega_m2*86400*x);
yfit_shear_new = shear_fit*sin(omega_m2*86400*x)+fittedmdl_shear.d;


figure(1);
clf;set(gcf,'Color','w')
hold on;
plot(x,y,'Color','k','LineWidth',1.5)
% l1=plot(fittedmdl_shear);
plot(xfit_shear,yfit_shear,'Color',yellow,'LineWidth',2);
% plot(xfit_shear,yfit_shear_new,'Color',red,'LineWidth',2);
grid on;grid minor
axis tight
set(gca,'FontSize',fontsize);
title('Linear-fit shear','Interpreter','latex')
ylabel('($s^{-1}$)','Interpreter','latex')
ylim([-2.5 2.5]*1e-3)
xlabel('Time (days)','Interpreter','latex')
legend('Linear-fit shear','Sinusoidal curve fitting','Interpreter','latex')
box on;
% xlim([10 14.2])
% xlim([0.5 3.7])
% print('-dpng','-r200',['figures/fig1' fname '.png']);





y = smooth_n2';
y=y(yidx);
mdl = fittype('a*sin(b*x+c)+d','indep','x','options',fo);
fittedmdl_N2 = fit(x,y,mdl,'start',[rand(),omega_m2*86400,pi/2,rand()])
xfit_N2 = x;
yfit_N2 = fittedmdl_N2(xfit_N2);


N2 = fittedmdl_N2.d;

alpha = fittedmdl_N2.a;

alpha_theory = N2*shear_fit/omega_m2*sind(topo)*cosd(topo);

alpha/alpha_theory

yfit_N2_new = alpha_theory*sin(omega_m2*86400*x+pi/2)+fittedmdl_N2.d;

Ri_curvefit = yfit_N2./(yfit_shear.^2)/(cosd(topo))^2;
Ri_curvefit_new = yfit_N2_new./(shear_theory.^2)/(cosd(topo))^2;




%%%%%% Compute predicted N^2, using the kinematic model




%%%%%%



figure(2);
clf;set(gcf,'Color','w')
hold on;
plot(x,n2(yidx),'Color',gray,'LineWidth',1.5)
plot(x,smooth_n2(yidx),'--','Color','k','LineWidth',1.5)
plot(xfit_N2,yfit_N2,'Color',yellow,'LineWidth',2);
% plot(xfit_N2,yfit_N2_new,'Color',red,'LineWidth',2);
grid on;grid minor
axis tight
set(gca,'FontSize',fontsize);
title('$\tilde N^2$','Interpreter','latex')
legend('Unsmoothed $\tilde N^2$','Smoothed $\tilde N^2$','Sinusoidal curve fitting','Interpreter','latex')
ylabel('($s^{-2}$)','Interpreter','latex')
ylim([-1 13]*1e-6)
xlabel('Time (days)','Interpreter','latex')
box on;
% xlim([10 14.2])
% xlim([0.5 3.7])
% print('-dpng','-r200',['figures/fig2' fname '.png']);

% plot(xfit_shear,yfit_shear/1e3*2+4e-6,'Color',red,'LineWidth',2);




Ri = n2./(shear_int.^2)/(cosd(topo)^2);
smooth_Ri = smooth_n2./(shear_int.^2)/(cosd(topo)^2);


dt = 10;
NTtide = 300;
omega = omega_m2;
Nt = NTtide/omega/dt;

shear_convec = cosd(topo)/sind(topo)*omega;
dt_ri = dt/1000;
tt_ri = dt_ri:dt_ri:Nt*dt;

shear = shear_fit;
Ri_inverse = (cosd(topo)*shear*cos(omega*tt_ri)).^2./(N2 - N2*sind(topo)*cosd(topo)/omega*shear*sin(omega*tt_ri));
Ri_min = 1/max(Ri_inverse) 


% clear Ri_inverse tt_ri

figure(3)
clf;set(gcf,'Color','w')
plot(x, 1./Ri(yidx),'Color',gray,'LineWidth',1.5);grid on;
hold on;
plot(x, 1./smooth_Ri(yidx),'--','Color','k','LineWidth',1.5);
plot(x, 1./Ri_min*ones(1,length(yidx)),'-','LineWidth',2);
% plot(x, 1./Ri_curvefit,'-','LineWidth',2);
% plot(x, 1./Ri_curvefit_new,'-','LineWidth',2);

set(gca,'FontSize',fontsize);
legend('Unsmoothed','Smoothed')
xlabel('Time (days)','Interpreter','latex')
grid on;grid minor
title('Inverse Richardson number','Interpreter','latex')
axis tight
% ylim([-0.5 5])
ylim([-0.25 4])
% xlim([10 14.2])
% xlim([0.5 3.7])
legend('Inverse Ri (unsmoothed $\tilde N^2$)', 'Inverse Ri (smoothed $\tilde N^2$)','Inverse minimum Ri, sinusoidal curve fitting','Interpreter','latex')

% print('-dpng','-r200',['figures/fig3' fname '.png']);


clear ustart uend fontsize time_n2 time_uw shear time_temp n2start n2idx uidx N2_zavg shear_linear smooth_N2_zavg



