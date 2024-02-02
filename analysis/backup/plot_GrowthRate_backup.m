

%%


clear;close all;
load('muk.mat')
plot_parm = Shear_parm;
NEXP = length(plot_parm);

 figure(1)
 clf;set(gcf,'Color','w','Position',[117 426 872 496])

for ne = 8:13
    load_all
    filename = [expdir expname '/RMSE_mean.mat'];
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
% plot_mitgcm_parm = [0.6:0.1:1.8]*1e-3;
plot_mitgcm_parm = [0.8 1.0 1.2 1.4 1.6 1.8]*1e-3;

gr_plot = gr_tt(8:13);
gr_plot = gr_plot.^(1/ncycle);

figure(2)
clf;
set(gcf,'Color','w')
p1 = plot(plot_parm,muk_mean_buoy,'LineWidth',2);
hold on;
p2 = plot(plot_parm,muk_mean_zeta,'LineWidth',2);
% p1 = plot(plot_parm,muk_mean_psi,'LineWidth',2);
plot(plot_parm,ones(1,NEXP),'k--','LineWidth',1)
% s1 = scatter(plot_mitgcm_parm(1:4),gr_plot(1:4),100,'Filled');
% s2 = scatter(plot_mitgcm_parm(5:8),gr_plot(5:8),100,'Filled');
% s3 = scatter(plot_mitgcm_parm(9:12),gr_plot(9:12),100,'Filled');
s3 = scatter(plot_mitgcm_parm,gr_plot,100,'Filled');

set(gca,'Fontsize',fontsize)
grid on;grid minor;
title('Mean Floquet exponents \mu (\lambda = 400 m)','Fontsize',fontsize+5)
xlabel('Velocity shear \Lambda (1/s)','Fontsize',fontsize+5)
legend([p1,p2,s3],'Floquet \mu_k^b', 'Floquet \mu_k^\zeta','\mu_k^b of MITgcm','Fontsize',fontsize+5,'Position',[0.2000 0.6655 0.1676 0.1747])
% legend([p1,s1],'Floquet \mu_k^\psi', '\mu_k^w of MITgcm','Fontsize',fontsize+5,'Position',[0.2000 0.6655 0.1676 0.1747])



