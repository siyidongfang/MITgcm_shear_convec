%%% 
%%% calc_GrowthRate_floquet.m
%%% 
%%% Calculate the instability growth rate of the Floquet analysis
%%% following Radko (2019), using kinetic energy and potential energy
%%% As a function of Richardson Number
clear all;close all;

N2const = 1e-6;
S_max = [0.6:0.2:1.8]*1e-3;
Ri_mean = N2const./(0.5*S_max.^2);
Ri_min = N2const./(S_max.^2);


% expdir = '../instability/exps_shear_L400_0.002/';
expdir = '../instability/exps_inviscid_noBC/';
% expdir = '../instability/exps_newbc/';

m1km =1000;
% Shear_parm = [0 1e-23 5e-5 0.0002:0.0001:0.0012 0.0014 0.0016 0.0018]; 
% Shear_parm = [0.6:0.2:2.2]*1e-3;
% % Shear_parm = [0.2 0.4 0.6 0.8 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.8 2.0 2.2 2.4 2.6]*1e-3; 
% Shear_parm = [0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4]*1e-3; 
Shear_parm = [0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4]*1e-3; 
NEXP = length(Shear_parm);
dz_group = 0.002*ones(1,NEXP);
dt_group = 0.002*ones(1,NEXP);

for ne = 1:NEXP
% for ne = [1 3]
    ne

    shear = Shear_parm(ne); 
    expname = ['H1500_topo4_Pt43200_N0.001_S' num2str(shear) '_lambda400_dz' num2str(dz_group(ne)) '_dt' num2str(dt_group(ne)) '/'];
    % expname = ['H1500_topo1e-100_Pt43200_N0.001_S' num2str(shear) '_lambda400_dz' num2str(dz_group(ne)) '_dt' num2str(dt_group(ne)) '_RK4+AB3'];

    load([expdir expname '/output.mat'],'re_buoy','re_zeta','re_psi','omega','dt','Ptide','Nt','Nr','zz','Nshear')
    % re_buoy = real(buoy);
    % re_zeta = real(zeta);
    % re_psi = real(psi);
    
    %%%%% Floquet stability
    oT = round(Ptide*omega/dt);% The time step after one tidal cycle

    tidx1 = 4*oT+1:5*oT-2;
    tidx2 = 5*oT+1:6*oT-2;
    
    zidx=2:Nshear;
    % zidx = 2:Nr-1;
    muk_psi = mean(abs(re_psi(tidx2,zidx)))./mean(abs(re_psi(tidx1,zidx)));
    muk_zeta = mean(abs(re_zeta(tidx2,zidx)))./mean(abs(re_zeta(tidx1,zidx)));
    muk_buoy = mean(abs(re_buoy(tidx2,zidx)))./mean(abs(re_buoy(tidx1,zidx)));

    % muk_radko = 

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

    % muk_max_radko(ne) = max(abs(muk_radko));
    % muk_mean_radko(ne) = mean(abs(muk_radko));
    % muk_rms_radko(ne) = rms(abs(muk_radko));

end

% log10_muk_buoy = log10(muk_mean_buoy);
% log10_muk_zeta = log10(muk_mean_zeta);
% log10_muk_psi = log10(muk_mean_psi);


save('muk_inviscid_noBC.mat','muk_mean_buoy','muk_mean_zeta','muk_mean_psi',...
               'muk_max_buoy','muk_max_zeta','muk_max_psi',...
               'muk_rms_buoy','muk_rms_zeta','muk_rms_psi',...
               'Shear_parm')



