%%%
%%% calc_FloquetMultipliers.m
%%%

clear;
expdir = 'exps_lambda_S0.0013_0.002/';

m1km =1000;

lambda_parm = [0.05 0.075 0.1 0.25 0.5 0.75 1 2.5 5 7.5 10 25 50 75 100 250 500 1000]*m1km;

% lambda_parm = [0.05 0.075 0.1 0.25 0.5 0.75 1 2.5 5 7.5 10 50 75 100]*m1km;
% lambda_parm = [0.25 0.5 0.75 1 1.2 1.5 1.7 2 2.5 3 3.5 4 5 7.5 10 25 50 75 100 250 500 1000 2000 5000 10000]*m1km;
kx_parm = 2*pi./lambda_parm;

NEXP = length(lambda_parm);

dz_group = 0.002*ones(1,NEXP);
dt_group = 0.002*ones(1,NEXP);

for ne = 1:NEXP
    lambda = lambda_parm(ne)
    expname = ['H1500_topo4_Pt43200_N0.001_S0.0013_lambda' num2str(lambda) '_dz' num2str(dz_group(ne)) '_dt' num2str(dt_group(ne)) '_RK4+AB3'];
    load([expdir expname '/output.mat'],'buoy','zeta','psi','omega','dt','Ptide','Nt','Nr','zz','Nshear')
    re_buoy = real(buoy);
    re_zeta = real(zeta);
    re_psi = real(psi);
    
    %%%%% Floquet stability
    oT = round(Ptide*omega/dt);% The time step after one tidal cycle
    tidx1 = 1*oT+1:4*oT-2;
    tidx2 = 3*oT+1:5*oT-2;
    zidx=2:Nshear;
    % zidx = 2:Nr-1;

    % tidx1 = 2*oT+1:2*oT+1;
    % tidx2 = 3*oT+1:3*oT+1;

    muk_psi = mean(abs(re_psi(tidx2,zidx)))./mean(abs(re_psi(tidx1,zidx)));
    muk_zeta = mean(abs(re_zeta(tidx2,zidx)))./mean(abs(re_zeta(tidx1,zidx)));
    muk_buoy = mean(abs(re_buoy(tidx2,zidx)))./mean(abs(re_buoy(tidx1,zidx)));
    % muk_psi = (abs(re_psi(tidx2,zidx)))./(abs(re_psi(tidx1,zidx)));
    % muk_zeta = (abs(re_zeta(tidx2,zidx)))./(abs(re_zeta(tidx1,zidx)));
    % muk_buoy = (abs(re_buoy(tidx2,zidx)))./(abs(re_buoy(tidx1,zidx)));


    muk_max_buoy(ne) = max(abs(muk_buoy));
    muk_mean_buoy(ne) = mean(abs(muk_buoy));
    muk_rms_buoy(ne) = rms(abs(muk_buoy));

    muk_max_zeta(ne) = max(abs(muk_zeta));
    muk_mean_zeta(ne) = mean(abs(muk_zeta));
    muk_rms_zeta(ne) = rms(abs(muk_zeta));

    muk_max_psi(ne) = max(abs(muk_psi));
    muk_mean_psi(ne) = mean(abs(muk_psi));
    muk_rms_psi(ne) = rms(abs(muk_psi));

end



%%
fontsize = 18;


% log10_muk_buoy = log10(muk_mean_buoy);
% log10_muk_zeta = log10(muk_mean_zeta);
% log10_muk_psi = log10(muk_mean_psi);

log10_muk_buoy = (muk_mean_buoy);
log10_muk_zeta = (muk_mean_zeta);
log10_muk_psi = (muk_mean_psi);

% log10_muk_buoy = (muk_max_buoy);
% log10_muk_zeta = (muk_max_zeta);
% log10_muk_psi = (muk_max_psi);

figure()
clf;
set(gcf,'Color','w','Position',[279 593 1072 397])
semilogx(lambda_parm,log10_muk_buoy,'LineWidth',2)
hold on;
semilogx(lambda_parm,log10_muk_zeta,'LineWidth',2)
semilogx(lambda_parm,log10_muk_psi,'LineWidth',2)
semilogx(lambda_parm,zeros(1,NEXP),'k--','LineWidth',1)
set(gca,'Fontsize',fontsize)
% ylim([-5 6])
grid on;grid minor;
title('Mean Floquet exponents log(\mu_k)','Fontsize',fontsize+5)
ylabel('log(\mu_k)','Fontsize',fontsize+5)
xlabel('Wavelength \lambda (m)','Fontsize',fontsize+5)
legend('\mu_k^b','\mu_k^\zeta','\mu_k^\psi')

figure()
clf;
set(gcf,'Color','w','Position',[279 593 1072 397])
semilogx(kx_parm,log10_muk_buoy,'LineWidth',2)
hold on;
semilogx(kx_parm,log10_muk_zeta,'LineWidth',2)
semilogx(kx_parm,log10_muk_psi,'LineWidth',2)
semilogx(kx_parm,zeros(1,NEXP),'k--','LineWidth',1)
set(gca,'Fontsize',fontsize)
% ylim([-5 6])
% xlim([1e-80 1])
grid on;grid minor;
title('Mean Floquet exponents log(\mu_k)','Fontsize',fontsize+5)
ylabel('log(\mu_k)','Fontsize',fontsize+5)
xlabel('Cross-isobath wavenumber k (m^{-1})','Fontsize',fontsize+5)
legend('\mu_k^b','\mu_k^\zeta','\mu_k^\psi')





