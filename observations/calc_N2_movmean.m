%%%%%%%% Compute the move-mean average of the sinusoidal fit shear
%%%%%%%% amplitude
clear all;
% load('MAVS2_N2.mat','smooth_N2_zavg','time_temp')
% time = 1:length(time_temp);N2 = smooth_N2_zavg; clear smooth_N2_zavg time_temp;
% isnanidx = 135753:length(N2);
% time = time(isnanidx);

load('MAVS1_N2.mat','smooth_N2_zavg','time_temp')
time = 1:length(time_temp);N2 = smooth_N2_zavg; clear smooth_N2_zavg time_temp;
isnanidx = 1:length(time);
time = time(isnanidx);

N2 = N2(isnanidx);
% %%% Interpolating NaN-s
% nanx = isnan(N2);
% t    = 1:numel(N2);
% N2(nanx) = interp1(t(~nanx), N2(~nanx), t(nanx));


fontsize = 16;
figure(1)
clf;set(gcf,'Color','w');
plot(time/86400,N2);
grid on;grid minor;
set(gca,'Fontsize',fontsize);
xlabel('Time (days)')
ylabel('N2 (1/s^2)')
title('N2 at MAVS1')
% ylim([-1.7 2.2]/1e3)
% xlim([0 90])


%%
%%% Curve fitting a sinusoidal function using a 4-tidal-cycle (~2 days) moving window
Nfit_cycle = 4;
Period_m2 = 12*3600+25*60;
omega_m2 = 2*pi/Period_m2;

x_all = time/86400;
y_all = N2;

Nend = length(N2)-round(Period_m2*Nfit_cycle)+1;
N2_fit = NaN*zeros(1,Nend);
N2_phase = NaN*zeros(1,Nend);
% N2_freq = NaN*zeros(1,Nend);

parfor n=1:Nend
% for n=2
    n
    Nstarty = n;
    Nendy = Nstarty + round(Period_m2*Nfit_cycle);
    if(Nendy>length(y_all))
        Nendy = length(y_all);
    end
    
    yidx = Nstarty:Nendy;
    y=y_all(yidx);
    x=x_all(yidx);
    % x=x-x(1);
    fo = fitoptions('Method','NonlinearLeastSquares');
    mdl = fittype('a*sin(12.1447*x+c)+d','indep','x','options',fo,coefficients=["a" "c" "d"]);
    fittedmdl_shear= fit(x',y,mdl,'start',[rand(),-0.44*pi+n/43200*2*pi,rand()]);
    % mdl = fittype('a*sin(b*x+c)+d','indep','x','options',fo,coefficients=["a" "b" "c" "d"]);
    % fittedmdl_shear= fit(x',y',mdl,'start',[rand(),omega_m2*86400,-0.44*pi+n*900/43200*2*pi,rand()]);
    

    xfit = x;
    yfit = fittedmdl_shear(xfit);
    
    N2_fit(n) = fittedmdl_shear.d;
    % N2_freq(n) = fittedmdl_shear.c;
    N2_phase(n) = fittedmdl_shear.c;
   
    % figure(2);clf;
    % set(gcf,'Color','w')
    % hold on;
    % plot(x,y,'LineWidth',1.5);
    % plot(xfit, yfit,'LineWidth',1.5);
    % grid on;grid minor;
    % xlabel('Time (days)');
    % set(gca,'FontSize',fontsize)
    % xlim([min(x) max(x)])
    % title('Observed N^2 (1/s^2)','FontSize',fontsize+5)

end



figure(3)
subplot(3,1,1)
plot(time(1:Nend)/86400,abs(N2_fit))
% subplot(3,1,2)
% plot(time(1:Nend),N2_freq)
subplot(3,1,3)
plot(time(1:Nend)/86400,N2_phase)


save('MAVS1_N2_movmean.mat')
