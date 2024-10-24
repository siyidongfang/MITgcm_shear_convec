%%%%%%%% Compute the move-mean average of the sinusoidal fit shear
%%%%%%%% amplitude

load('MAVS2_LinearShear.mat','time','shear_linear','topo')
isnanidx = 32:8574;
% load('MAVS1_LinearShear.mat','time','shear_linear','topo')
% isnanidx = 40:8756;

%%% Interpolating NaN-s
shear_linear = shear_linear(isnanidx);
time = time(isnanidx);
nanx = isnan(shear_linear);
t    = 1:numel(shear_linear);
shear_linear(nanx) = interp1(t(~nanx), shear_linear(~nanx), t(nanx));


fontsize = 16;
figure(1)
clf;set(gcf,'Color','w');
plot(time/86400,shear_linear);
grid on;grid minor;
set(gca,'Fontsize',fontsize);
xlabel('Time (days)')
ylabel('shear (1/s)')
title('Linear-fit shear at MAVS2')
ylim([-1.7 2.2]/1e3)
xlim([0 90])


%%
%%% Curve fitting a sinusoidal function using a 4-tidal-cycle (~2 days) moving window
Nfit_cycle = 4;
Period_m2 = 12*3600+25*60;
omega_m2 = 2*pi/Period_m2;

x_all = time/86400;
y_all = shear_linear;

Nend = length(shear_linear)-round(Period_m2/3600*4*Nfit_cycle)+1;
shear_fit = NaN*zeros(1,Nend);
shear_phase = NaN*zeros(1,Nend);
shear_freq = NaN*zeros(1,Nend);

parfor n=1:Nend
    Nstarty = n;
    Nendy = Nstarty + round(Period_m2/3600*4*Nfit_cycle);
    if(Nendy>length(y_all))
        Nendy = length(y_all);
    end
    
    yidx = Nstarty:Nendy;
    y=y_all(yidx);
    x=x_all(yidx);
    % x=x-x(1);
    fo = fitoptions('Method','NonlinearLeastSquares');
    mdl = fittype('a*sin(12.1447*x+c)+d','indep','x','options',fo,coefficients=["a" "c" "d"]);
    fittedmdl_shear= fit(x',y',mdl,'start',[rand(),-0.44*pi+n*900/43200*2*pi,rand()]);
    % mdl = fittype('a*sin(b*x+c)+d','indep','x','options',fo,coefficients=["a" "b" "c" "d"]);
    % fittedmdl_shear= fit(x',y',mdl,'start',[rand(),omega_m2*86400,-0.44*pi+n*900/43200*2*pi,rand()]);
    

    if(abs(fittedmdl_shear.a)<0.2e-3)
        n
        fittedmdl_shear1= fit(x',y',mdl,'start',[rand(),0*pi,rand()]);
        fittedmdl_shear2= fit(x',y',mdl,'start',[rand(),0.25*pi,rand()]);
        fittedmdl_shear3= fit(x',y',mdl,'start',[rand(),0.5*pi,rand()]);
        fittedmdl_shear4= fit(x',y',mdl,'start',[rand(),0.75*pi,rand()]);
        fittedmdl_shear5= fit(x',y',mdl,'start',[rand(),1*pi,rand()]);
        fittedmdl_shear6= fit(x',y',mdl,'start',[rand(),1.25*pi,rand()]);
        fittedmdl_shear7= fit(x',y',mdl,'start',[rand(),1.75*pi,rand()]);
        fittedmdl_shear8= fit(x',y',mdl,'start',[rand(),1.75*pi,rand()]);

        % fittedmdl_shear1= fit(x',y',mdl,'start',[rand(),omega_m2*86400,0*pi,rand()]);
        % fittedmdl_shear2= fit(x',y',mdl,'start',[rand(),omega_m2*86400,0.25*pi,rand()]);
        % fittedmdl_shear3= fit(x',y',mdl,'start',[rand(),omega_m2*86400,0.5*pi,rand()]);
        % fittedmdl_shear4= fit(x',y',mdl,'start',[rand(),omega_m2*86400,0.75*pi,rand()]);
        % fittedmdl_shear5= fit(x',y',mdl,'start',[rand(),omega_m2*86400,1*pi,rand()]);
        % fittedmdl_shear6= fit(x',y',mdl,'start',[rand(),omega_m2*86400,1.25*pi,rand()]);
        % fittedmdl_shear7= fit(x',y',mdl,'start',[rand(),omega_m2*86400,1.75*pi,rand()]);
        % fittedmdl_shear8= fit(x',y',mdl,'start',[rand(),omega_m2*86400,1.75*pi,rand()]);

        a1 = fittedmdl_shear1.a;
        a2 = fittedmdl_shear2.a;
        a3 = fittedmdl_shear3.a;
        a4 = fittedmdl_shear4.a;
        a5 = fittedmdl_shear5.a;
        a6 = fittedmdl_shear6.a;
        a7 = fittedmdl_shear7.a;
        a8 = fittedmdl_shear8.a;

        [maxa idxa] = max(abs([a1 a2 a3 a4 a5 a6 a7 a8]));
        fittedmdl_shear = fit(x',y',mdl,'start',[rand(),(idxa-1)*0.25*pi,rand()]);
        % fittedmdl_shear = fit(x',y',mdl,'start',[rand(),omega_m2*86400,(idxa-1)*0.25*pi,rand()]);
    end

    
    xfit = x;
    yfit = fittedmdl_shear(xfit);
    
    shear_fit(n) = fittedmdl_shear.a;
    shear_phase(n) = fittedmdl_shear.c;
    % shear_freq(n) =  fittedmdl_shear.b;
    % shear_theory = shear_fit*sin(omega_m2*86400*x+shear_phase);
    % yfit_shear_new = shear_fit*sin(omega_m2*86400*x+shear_phase)+fittedmdl_shear.d;
    
    % figure(2);clf;
    % set(gcf,'Color','w')
    % hold on;
    % % plot(xfit, shear_theory)
    % % plot(xfit, yfit_shear_new)
    % plot(x,y,'LineWidth',1.5);
    % plot(xfit, yfit,'LineWidth',1.5);
    % grid on;grid minor;
    % xlabel('Time (days)');
    % set(gca,'FontSize',fontsize)
    % xlim([min(x) max(x)])
    % title('Observed velocity shear (1/s)','FontSize',fontsize+5)

end



figure(3)
subplot(3,1,1)
plot(time(1:Nend)/86400,abs(shear_fit))
% subplot(3,1,2)
% plot(time(1:Nend),shear_freq)
subplot(3,1,3)
plot(time(1:Nend)/86400,shear_phase)


save('MAVS2_shear_movmean.mat')
