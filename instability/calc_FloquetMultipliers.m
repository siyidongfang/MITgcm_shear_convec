%%%
%%% calc_FloquetMultipliers.m
%%%

clear;

m1km =1000;
lambda_parm = [0.005 0.0075 0.01 0.025 0.05 0.075 0.1 0.25 0.5 0.75 1 2.5 5 7.5 10 25 50 75 100 250 500 1000]*m1km;
dz_group = 0.01*ones(1,22);
dt_group = 0.01*ones(1,22);

dt_group(4) = 0.003;
dt_group(3) = 0.002;
dt_group(2) = 0.002;
dt_group(1) = 0.001;


for ne = 1:22
    lambda = lambda_parm(ne)
    expname = ['topo4_Pt43200_N0.001_S0.001_L' num2str(lambda) '_dz' num2str(dz_group(ne)) '_dt' num2str(dt_group(ne))];
    
    load(['exps/' expname '/output.mat'],'buoy','zeta','psi','omega','dt','Ptide','Nt','Nr','zz','Nshear')
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
set(gcf,'Color','w')
semilogx(lambda_parm,ploty,'LineWidth',2)
hold on;
semilogx(lambda_parm,zeros(1,22),'k--','LineWidth',1)
set(gca,'Fontsize',fontsize)
% ylim([0 60])
grid on;grid minor;
title('Mean Floquet exponents log(\mu_k)','Fontsize',fontsize+5)
ylabel('log(\mu_k)','Fontsize',fontsize+5)
xlabel('\lambda (m)','Fontsize',fontsize+5)





