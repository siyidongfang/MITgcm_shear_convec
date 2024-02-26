
clear;close all;
load('exps_NoStressBottom0.01.mat')
% load('muk.mat')
plot_parm = Shear_parm;
NEXP = length(plot_parm);

Pr = 1;

addpath    /Users/ysi/MITgcm_shear_convec/analysis

for ne = 1:7
    load_all
    filename = [expdir expname '/RMSE_mean.mat'];
    load(filename);

    %%% Calculate the growth rate (the second tidal cycle divided by the first)
    tidx2 = 12+1:24;
    tidx1 = 0+1:12;
    % tidx2 =25:36;
    % tidx1 =12+1:24;
    % tidx2 =37:48;
    % tidx1 =25:36;

    ncycle = (tidx2(1)-tidx1(1))/12;
    
    % gr_tt(ne) = mean(div_tt_zavg(tidx2)./div_tt_zavg(tidx1));
    % gr_uu(ne) = mean(div_uu_zavg(tidx2)./div_uu_zavg(tidx1));
    % gr_ww(ne) = mean(div_ww_zavg(tidx2)./div_ww_zavg(tidx1));
    % 
    gr_tt(ne) = max(div_tt_zavg(tidx2)./div_tt_zavg(tidx1));
    gr_uu(ne) = max(div_uu_zavg(tidx2)./div_uu_zavg(tidx1));
    gr_ww(ne) = max(div_ww_zavg(tidx2)./div_ww_zavg(tidx1));

   
    % gr_radko(ne) = mean(mean((div_tt(tidx2,:).^2)./(div_tt(tidx1,:).^2)*Pr/4 ...
    %     + (div_uu(tidx2,:).^2)./(div_uu(tidx1,:).^2)/4 ...
    %     + (div_ww(tidx2,:).^2)./(div_ww(tidx1,:).^2)/4)); ...



end

% S_max_gcm = [0.8:0.1:1.9]*1e-3;

S_max_gcm = [0.4:0.2:1.6]*1e-3;

gr_plot = gr_tt;
gr_plot = gr_plot.^(1/ncycle);


%%
figure(1)
clf;
set(gcf,'Color','w')
p1 = plot(plot_parm,muk_max_buoy,'LineWidth',2);
hold on;
p2 = plot(plot_parm,muk_max_zeta,'LineWidth',2);
plot(plot_parm,ones(1,NEXP),'k--','LineWidth',1)
s3 = scatter(S_max_gcm,gr_plot,100,'Filled');

set(gca,'Fontsize',fontsize)
grid on;grid minor;
title('Max Floquet exponents \mu (\lambda = 400 m)','Fontsize',fontsize+5)
xlabel('Velocity shear \Lambda (1/s)','Fontsize',fontsize+5)
legend([p1,p2,s3],'Floquet \mu_k^b', 'Floquet \mu_k^\zeta','\mu_k^b of MITgcm','Fontsize',fontsize+5,'Position',[0.2000 0.6655 0.1676 0.1747])


%%
%%% Calculate the growth rate defined in Radko 2019 as a function of mean
%%% Richardson Number
N2const = 1e-6;
Ri_mean_gcm = N2const./(0.5*S_max_gcm.^2);
Ri_mean_flo = N2const./(0.5*Shear_parm.^2);

lam_gcm = log(gr_plot)/2;
lam_flo_mean = log(muk_mean_buoy)/2;
lam_flo_max  = log(muk_max_buoy)/2;


cc1 = 0.3553;
cc2 = -0.01467;

Ri_radko = 0.05:0.05:15;
lam_norm_radko = cc1./sqrt(Ri_radko)+cc2./Ri_radko; 



figure(2)
clf;set(gcf,'Color','w')
Lflo_max = plot(Ri_mean_flo,lam_flo_max,'LineWidth',2,'Color',blue);
hold on;
Lflo_mean = plot(Ri_mean_flo,lam_flo_mean,':','LineWidth',2,'Color',blue);
Lradko = plot(Ri_radko,lam_norm_radko,'--','LineWidth',2);
Sgcm = scatter(Ri_mean_gcm,lam_gcm,100,'Filled');

set(gca,'Fontsize',fontsize)
grid on;grid minor;
xlabel('$\overline{R_i}$','Interpreter','latex','Fontsize',fontsize+5)
lg1 = legend([Lflo_max,Lflo_mean,Sgcm,Lradko],...
    'Floquet (max)','Floquet (mean)','MITgcm', 'Radko (2019)',...
    'Fontsize',fontsize+3,'Position', [0.3789 0.6219 0.3170 0.2048]);
title('Growth Rate \lambda')












