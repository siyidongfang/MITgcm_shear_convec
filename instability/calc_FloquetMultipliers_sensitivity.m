%%%
%%% calc_FloquetMultipliers.m
%%%

clear;
expdir = 'exps_shear_L400_0.002/';

m1km =1000;

topo_parm = [1e-20 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
% Shear_parm = [0 1e-20 0.05 0.1:0.1:1 1.2 1.4 1.6 1.8 2 2.2 2.4]*1e-3; 
% Shear_parm = [0 1e-20 0.05 0.1:0.1:1]*1e-3; 
Shear_parm = [0 1e-23 5e-5 0.0005 0.0006 0.0007 0.0008 0.0009 0.001 0.0011 0.0012 0.0014 0.0016 0.0018 0.002]; 
N_parm = [1e-20 0.01 0.05 0.1 0.25 0.5 0.75 1 2 3 4]*1e-3;
dz_parm = [0.1:0.1:2]*0.01;
dt_parm = [0.001 0.002 0.005 0.007 0.01];


plot_parm = Shear_parm

NEXP = length(plot_parm);

dz_group = 0.002*ones(1,NEXP);
dt_group = 0.002*ones(1,NEXP);

for ne = 1:NEXP
    ne
    % topo = topo_parm(ne)
    % expname = ['H1500_topo' num2str(topo) '_Pt43200_N0.001_S0.001_lambda1000_dz' num2str(dz_group(ne)) '_dt' num2str(dt_group(ne))];
    
    shear = Shear_parm(ne); 
    expname = ['H1500_topo4_Pt43200_N0.001_S' num2str(shear) '_lambda400_dz' num2str(dz_group(ne)) '_dt' num2str(dt_group(ne)) '_RK4+AB3'];

    % N = N_parm(ne); 
    % expname = ['H1500_topo4_Pt43200_N' num2str(N) '_S0.001_lambda1000_dz' num2str(dz_group(ne)) '_dt' num2str(dt_group(ne))];

    % dz = 0.002;
    % dt = dt_parm(ne)
    % dz = dz_parm(ne)
    % dt = dz;
    % expname = ['H1500_topo4_Pt43200_N0.001_S0.001_lambda400_dz' num2str(dz) '_dt' num2str(dt) '_RK4+AB3'];


    load([expdir expname '/output.mat'],'buoy','zeta','psi','omega','dt','Ptide','Nt','Nr','zz','Nshear')
    re_buoy = real(buoy);
    re_zeta = real(zeta);
    re_psi = real(psi);
    
    %%%%% Floquet stability
    oT = round(Ptide*omega/dt);% The time step after one tidal cycle
    % tidx1 = 12*oT+1:16*oT;
    % tidx2 = 16*oT+1:20*oT;
    tidx1 = 2*oT+1:3*oT-2;
    tidx2 = 3*oT+1:4*oT-2;
    
    % zidx=2:Nshear;
    zidx = 2:Nr-1;
    muk_psi = mean(abs(re_psi(tidx2,zidx)))./mean(abs(re_psi(tidx1,zidx)));
    muk_zeta = mean(abs(re_zeta(tidx2,zidx)))./mean(abs(re_zeta(tidx1,zidx)));
    muk_buoy = mean(abs(re_buoy(tidx2,zidx)))./mean(abs(re_buoy(tidx1,zidx)));
    % muk_psi = mean(re_psi(tidx2,zidx))./mean(re_psi(tidx1,zidx));
    % muk_zeta = mean(re_zeta(tidx2,zidx))./mean(re_zeta(tidx1,zidx));
    % muk_buoy = mean(re_buoy(tidx2,zidx))./mean(re_buoy(tidx1,zidx));
    
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


log10_muk_buoy = log10(muk_mean_buoy);
log10_muk_zeta = log10(muk_mean_zeta);
log10_muk_psi = log10(muk_mean_psi);

% log10_muk_buoy = (muk_mean_buoy);
% log10_muk_zeta = (muk_mean_zeta);
% log10_muk_psi = (muk_mean_psi);

figure()
clf;
set(gcf,'Color','w')
plot(plot_parm,log10_muk_buoy,'LineWidth',2)
hold on;
plot(plot_parm,log10_muk_zeta,'LineWidth',2)
plot(plot_parm,log10_muk_psi,'LineWidth',2)
plot(plot_parm,zeros(1,NEXP),'k--','LineWidth',1)
set(gca,'Fontsize',fontsize)
grid on;grid minor;
title('Mean Floquet exponents log(\mu) (\lambda = 400 m)','Fontsize',fontsize+5)
% ylabel('log(\mu)','Fontsize',fontsize+5)
% xlabel('topographic slope \theta','Fontsize',fontsize+5)
xlabel('Velocity shear \Lambda (1/s)','Fontsize',fontsize+5)
% xlabel('Buoyancy frequency N (1/s)','Fontsize',fontsize+5)
% xlabel('dt','Fontsize',fontsize+5)

legend('\mu_k^b','\mu_k^\zeta','\mu_k^\psi')

