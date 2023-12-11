%%%
%%% calc_FloquetMultipliers.m
%%%

clear;
expdir = 'lambda200_shear/';

m1km =1000;

topo_parm = [1e-20 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
% Shear_parm = [1e-20 0.05 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2 2.2 2.4 2.6]*1e-3; %%% lambda = 1000
Shear_parm = [1e-20 0.05 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2 2.2 2.4]*1e-3; %%% lambda = 1000
N_parm = [1e-20 0.01 0.05 0.1 0.25 0.5 0.75 1 2 3 4]*1e-3;

NEXP = length(Shear_parm);

dz_group = 0.01*ones(1,NEXP);
dt_group = 0.01*ones(1,NEXP);
% dt_group(9) = 0.002;

for ne = 1:NEXP
    % topo = topo_parm(ne)
    % expname = ['H1500_topo' num2str(topo) '_Pt43200_N0.001_S0.001_lambda1000_dz' num2str(dz_group(ne)) '_dt' num2str(dt_group(ne))];
    
    shear = Shear_parm(ne); 
    expname = ['H1500_topo4_Pt43200_N0.001_S' num2str(shear) '_lambda200_dz' num2str(dz_group(ne)) '_dt' num2str(dt_group(ne))];

    % N = N_parm(ne); 
    % expname = ['H1500_topo4_Pt43200_N' num2str(N) '_S0.001_lambda1000_dz' num2str(dz_group(ne)) '_dt' num2str(dt_group(ne))];

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
    
    % sum(abs(muk_psi)>1)
    % sum(abs(muk_zeta)>1)
    % sum(abs(muk_buoy)>1)
    % figure(1);
    % clf;
    % plot(muk_psi,zz(zidx))
    % hold on;
    % plot(muk_zeta,zz(zidx))
    % plot(muk_buoy,zz(zidx))
    % legend('psi','zeta','buoy')

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
set(gcf,'Color','w')
plot(Shear_parm,log10_muk_buoy,'LineWidth',2)
hold on;
plot(Shear_parm,log10_muk_zeta,'LineWidth',2)
plot(Shear_parm,log10_muk_psi,'LineWidth',2)
% plot(topo_parm,zeros(1,NEXP),'k--','LineWidth',1)
set(gca,'Fontsize',fontsize)
grid on;grid minor;
title('Maximum Floquet exponents log(\mu) (\lambda = 200 m)','Fontsize',fontsize+5)
ylabel('log(\mu)','Fontsize',fontsize+5)
% xlabel('topographic slope \theta','Fontsize',fontsize+5)
xlabel('Velocity shear \Lambda (1/s)','Fontsize',fontsize+5)
% xlabel('Buoyancy frequency N (1/s)','Fontsize',fontsize+5)
legend('\mu_k^b','\mu_k^\zeta','\mu_k^\psi')

