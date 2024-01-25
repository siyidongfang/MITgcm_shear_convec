

%% Floquet



expdir = '../instability/exps_shear_L400_0.002/';

m1km =1000;
Shear_parm = [0 1e-23 5e-5 0.0002:0.0001:0.0012 0.0014 0.0016 0.0018]; 
NEXP = length(plot_parm);
dz_group = 0.002*ones(1,NEXP);
dt_group = 0.002*ones(1,NEXP);

for ne = 1:NEXP
    ne

    shear = Shear_parm(ne); 
    expname = ['H1500_topo4_Pt43200_N0.001_S' num2str(shear) '_lambda400_dz' num2str(dz_group(ne)) '_dt' num2str(dt_group(ne)) '_RK4+AB3'];

    load([expdir expname '/output.mat'],'buoy','zeta','psi','omega','dt','Ptide','Nt','Nr','zz','Nshear')
    re_buoy = real(buoy);
    re_zeta = real(zeta);
    re_psi = real(psi);
    
    %%%%% Floquet stability
    oT = round(Ptide*omega/dt);% The time step after one tidal cycle
    % tidx1 = 12*oT+1:16*oT;
    % tidx2 = 16*oT+1:20*oT;

    tidx1 = 3*oT+1:4*oT-2;
    tidx2 = 4*oT+1:5*oT-2;
    % tidx1 = 2*oT+1:3*oT;
    % tidx2 = 3*oT+1:4*oT;
    % tidx1 = 0*oT+1:1*oT;
    % tidx2 = 1*oT+1:2*oT;
    
    zidx=2:Nshear;
    % zidx = 2:Nr-1;
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

% log10_muk_buoy = log10(muk_mean_buoy);
% log10_muk_zeta = log10(muk_mean_zeta);
% log10_muk_psi = log10(muk_mean_psi);


save('muk.mat','muk_mean_buoy','muk_mean_zeta','muk_mean_psi',...
               'muk_max_buoy','muk_max_zeta','muk_max_psi',...
               'muk_rms_buoy','muk_rms_zeta','muk_rms_psi',...
   'Shear_parm')


%%


clear;close all;
load('muk.mat')
plot_parm = Shear_parm;
NEXP = length(plot_parm);

 figure(1)
 clf;set(gcf,'Color','w','Position',[117 426 872 496])

for ne = 1:13
    load_all
    filename = [expdir expname '/RMSE_inst.mat'];
    load(filename);

    %%% Calculate the growth rate (the second tidal cycle divided by the first)
    tidx2 = 12+1:24;
    tidx1 = 0+1:12;

    % tidx2 = 36+1:48;
    % tidx1 = 24+1:36;

    % tidx2 = 24+1:36;
    % tidx1 = 12+1:24;

    % tidx2 = 12+1:24+6;
    % tidx1 = 0+1:12+6;

    % tidx2 = 48+1:60;
    % tidx1 = 36+1:48;

    % tidx2 = 48+1:60;
    % tidx1 = 0+1:12;  

    % tidx2 = 36+1:48;
    % tidx1 = 12+1:24;

    ncycle = (tidx2(1)-tidx1(1))/12;


    gr_tt(ne) = mean(div_tt_zavg(tidx2)./div_tt_zavg(tidx1));
    gr_uu(ne) = mean(div_uu_zavg(tidx2)./div_uu_zavg(tidx1));
    % gr_vv(ne) = mean(div_vv_zavg(tidx2)./div_vv_zavg(tidx1));
    gr_ww(ne) = mean(div_ww_zavg(tidx2)./div_ww_zavg(tidx1));
   
    % figure(1)
    % hold on;
    % plot(time_h,div_tt_zavg,'LineWidth',2);
end



    % set(gca,'Fontsize',fontsize)
    % xlabel('Time (hours)')
    % title('Temperature RMSE averaged over the bottom shear layer')
    % % title('Velocity u RMSE averaged over the bottom shear layer')
    % ylabel('(degC)')
    % % ylabel('(m/s)')
    % grid on;grid minor;
    % legend(...
    %     '\Lambda = 0.6\times10^{-3} s^{-1}',...
    %     '\Lambda = 0.7\times10^{-3} s^{-1}',...
    %     '\Lambda = 0.8\times10^{-3} s^{-1}',...
    %     '\Lambda = 0.9\times10^{-3} s^{-1}',...
    %     '\Lambda = 1.0\times10^{-3} s^{-1}',...
    %     '\Lambda = 1.1\times10^{-3} s^{-1}',...
    %     '\Lambda = 1.2\times10^{-3} s^{-1}',...
    %     '\Lambda = 1.3\times10^{-3} s^{-1}',...
    %     '\Lambda = 1.4\times10^{-3} s^{-1}',...
    %     '\Lambda = 1.5\times10^{-3} s^{-1}',...
    %     '\Lambda = 1.7\times10^{-3} s^{-1}',...
    %     '\Lambda = 1.8\times10^{-3} s^{-1}',...
    %     'Fontsize',fontsize,'Position', [0.1711 0.6063 0.2024 0.2288])


    %%
plot_mitgcm_parm = [0.6:0.1:1.8]*1e-3;

gr_plot = gr_tt;
gr_plot = gr_plot.^(1/ncycle);

figure(2)
clf;
set(gcf,'Color','w')
% p1 = plot(plot_parm,muk_mean_buoy,'LineWidth',2);
hold on;
% p1 = plot(plot_parm,muk_mean_zeta,'LineWidth',2);
p1 = plot(plot_parm,muk_mean_psi,'LineWidth',2);
plot(plot_parm,ones(1,NEXP),'k--','LineWidth',1)
s1 = scatter(plot_mitgcm_parm(1:4),gr_plot(1:4),100,'Filled');
s2 = scatter(plot_mitgcm_parm(5:8),gr_plot(5:8),100,'Filled');
s3 = scatter(plot_mitgcm_parm(9:12),gr_plot(9:12),100,'Filled');

set(gca,'Fontsize',fontsize)
grid on;grid minor;
title('Mean Floquet exponents \mu (\lambda = 400 m)','Fontsize',fontsize+5)
xlabel('Velocity shear \Lambda (1/s)','Fontsize',fontsize+5)
legend([p1,s1],'Floquet \mu_k^b', '\mu_k^b of MITgcm','Fontsize',fontsize+5,'Position',[0.2000 0.6655 0.1676 0.1747])
% legend([p1,s1],'Floquet \mu_k^\psi', '\mu_k^w of MITgcm','Fontsize',fontsize+5,'Position',[0.2000 0.6655 0.1676 0.1747])



