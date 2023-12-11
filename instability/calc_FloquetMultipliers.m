%%%
%%% calc_FloquetMultipliers.m
%%%

clear;
expdir = 'h500_dpsidz0/';

m1km =1000;
% lambda_parm = [0.001 0.0025 0.005 0.0075 0.01 0.025 0.05 0.075 0.1 0.25 0.5 0.75 1 2.5 5 7.5 10 25 50 75 100 250 500 1000]*m1km;
% lambda_parm = [0.0025 0.005 0.0075 0.01 0.025 0.05 0.075 0.1 0.25 0.5 0.75 1 2.5 5 7.5 10 25 50 75 100 250 500 1000 1e17 1.0e37 1e57 1e77]*m1km;
% lambda_parm = [0.05 0.075 0.1 0.25 0.5 0.75 1 2.5 5 7.5 10 25 50 75 100 250 500 1000]*m1km;
% lambda_parm = [0.05 0.075 0.1 0.25 0.5 0.75 1 1.2 1.5 1.7 2 2.5 3 3.5 4 5 7.5 10 25 50 75 100 250 500 1000 0]*m1km;
% lambda_parm = [0.1 0.25 0.5 0.75 1 2.5 5 7.5 10 25]*m1km;
lambda_parm = [250 500 750 1000 2000 2500 5000]
% lambda_parm = [0.001 0.0025 0.005 0.0075 0.01 0.025 0.05 0.075 0.1 0.25 0.5 0.75 1 1.2 1.5 1.7 2 2.5 3 3.5 4 5 7.5 10 25 50 75 100 250 500 1000 2000 5000 10000]*m1km;
% lambda_parm = [0.001 0.0025 0.005 0.0075 0.01 0.025 0.05 0.075 0.1 0.25 0.5 0.75 1 2.5 5 7.5 10 25 50 75 100 250 500 1000]*m1km;

kx_parm = 2*pi./lambda_parm;

NEXP = length(lambda_parm);

dz_group = 0.01*ones(1,NEXP);
dt_group = 0.01*ones(1,NEXP);

% dt_group(1) = 0.0002;
% dt_group(2) = 0.0003;
% dt_group(3) = 0.001;
% 
% dt_group(4:6) = 0.001;
% dt_group(7) = 0.005;




for ne = 1:NEXP
    lambda = lambda_parm(ne)
    expname = ['H500_topo4_Pt43200_N0.001_S0.001_lambda' num2str(lambda) '_dz' num2str(dz_group(ne)) '_dt' num2str(dt_group(ne))];
    % expname = ['topo4_Pt43200_N0.001_S0.001_L' num2str(lambda) '_dz' num2str(dz_group(ne)) '_dt' num2str(dt_group(ne))];
    load([expdir expname '/output.mat'],'buoy','zeta','psi','omega','dt','Ptide','Nt','Nr','zz','Nshear')
    re_buoy = real(buoy);
    re_zeta = real(zeta);
    re_psi = real(psi);
    
    %%%%% Floquet stability
    oT = round(Ptide*omega/dt);% The time step after one tidal cycle
    tidx = 100:100;
    zidx=2:Nshear;
    % zidx = 2:Nr-1;
    muk_psi = mean(re_psi(oT+tidx,zidx))./mean(re_psi(tidx,zidx));
    muk_zeta = mean(re_zeta(oT+tidx,zidx))./mean(re_zeta(tidx,zidx));
    muk_buoy = mean(re_buoy(oT+tidx,zidx))./mean(re_buoy(tidx,zidx));
    
    % sum(muk_psi>1)
    % sum(muk_zeta>1)
    % sum(muk_buoy>1)
    % figure(1);
    % clf;
    % plot(muk_psi,zz(zidx))
    % hold on;
    % plot(muk_zeta,zz(zidx))
    % plot(muk_buoy,zz(zidx))
    % legend('psi','zeta','buoy')
    
    % muk_max(ne) = max(abs([muk_buoy muk_zeta muk_psi]));
    % muk_mean(ne) = mean(abs([muk_buoy muk_zeta muk_psi]));
    % muk_rms(ne) = rms(abs([muk_buoy muk_zeta muk_psi]));

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


log10_muk_buoy = log10(muk_max_buoy);
log10_muk_zeta = log10(muk_max_zeta);
log10_muk_psi = log10(muk_max_psi);

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





