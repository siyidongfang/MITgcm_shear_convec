%%%
%%% calc_FloquetMultipliers.m
%%%

clear;
expdir = 'exps/';

m1km =1000;
% lambda_parm = [0.001 0.0025 0.005 0.0075 0.01 0.025 0.05 0.075 0.1 0.25 0.5 0.75 1 2.5 5 7.5 10 25 50 75 100 250 500 1000]*m1km;
% lambda_parm = [0.0025 0.005 0.0075 0.01 0.025 0.05 0.075 0.1 0.25 0.5 0.75 1 2.5 5 7.5 10 25 50 75 100 250 500 1000 1e17 1.0e37 1e57 1e77]*m1km;
lambda_parm = [0.05 0.075 0.1 0.25 0.5 0.75 1 2.5 5 7.5 10 25 50 75 100 250 500 1000 1e10 1e50 1e100]*m1km;
kx_parm = 2*pi./lambda_parm;

NEXP = length(lambda_parm);

dz_group = 0.01*ones(1,NEXP);
dt_group = 0.01*ones(1,NEXP);

dt_group(1) = 0.005;
dz_group(1) = 0.007;
% n=1;
% dt_group(6-n) = 0.003;
% dt_group(5-n) = 0.002;
% dt_group(4-n) = 0.002;
% dt_group(3-n) = 0.001;
% dt_group(2-n) = 0.0003;
% % dt_group(1) = 0.0002;

for ne = 1:NEXP
    lambda = lambda_parm(ne)
    if(ne==25)
        lambda = 1e40;
    end
    expname = ['H1500_topo4_Pt43200_N0.001_S0.001_lambda' num2str(lambda) '_dz' num2str(dz_group(ne)) '_dt' num2str(dt_group(ne))];
    
    load([expdir expname '/output.mat'],'buoy','zeta','psi','omega','dt','Ptide','Nt','Nr','zz','Nshear')
    re_buoy = real(buoy);
    re_zeta = real(zeta);
    re_psi = real(psi);
    
    %%%%% Floquet stability
    oT = round(Ptide*omega/dt);% The time step after one tidal cycle
    ntimestep = 100;
    zidx=2:Nshear;
    % zidx = 2:Nr-1;
    muk_psi = re_psi(oT+ntimestep,zidx)./re_psi(ntimestep,zidx);
    muk_zeta = re_zeta(oT+ntimestep,zidx)./re_zeta(ntimestep,zidx);
    muk_buoy = re_buoy(oT+ntimestep,zidx)./re_buoy(ntimestep,zidx);
    
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
    
    muk_max(ne) = max(abs([muk_buoy muk_zeta muk_psi]));
    muk_mean(ne) = mean(abs([muk_buoy muk_zeta muk_psi]));
    muk_rms(ne) = mean(abs([muk_buoy muk_zeta muk_psi]));

end




%%
fontsize = 18;


ploty = log10(muk_mean);

figure(2)
% clf;
set(gcf,'Color','w')
semilogx(lambda_parm,ploty,'LineWidth',2)
hold on;
semilogx(lambda_parm,zeros(1,NEXP),'k--','LineWidth',1)
set(gca,'Fontsize',fontsize)
% ylim([-5 6])
grid on;grid minor;
title('Mean Floquet exponents log(\mu_k)','Fontsize',fontsize+5)
ylabel('log(\mu_k)','Fontsize',fontsize+5)
xlabel('Wavelength \lambda (m)','Fontsize',fontsize+5)


figure(3)
% clf;
set(gcf,'Color','w')
semilogx(kx_parm,ploty,'LineWidth',2)
hold on;
semilogx(kx_parm,zeros(1,NEXP),'k--','LineWidth',1)
set(gca,'Fontsize',fontsize)
% ylim([-5 6])
xlim([1e-80 1])
grid on;grid minor;
title('Mean Floquet exponents log(\mu_k)','Fontsize',fontsize+5)
ylabel('log(\mu_k)','Fontsize',fontsize+5)
xlabel('Cross-isobath wavenumber k (m^{-1})','Fontsize',fontsize+5)





