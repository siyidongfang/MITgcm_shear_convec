

clear;close all;

 figure(1)
 clf;set(gcf,'Color','w','Position',[117 426 872 496])

for ne = 1:9
    load_all
    filename = [expdir expname '/RMSE_inst.mat'];
    load(filename);

    %%% Calculate the growth rate (the second tidal cycle divided by the first)
    % tidx2 = 12+1:24;
    % tidx1 = 0+1:12;

    tidx2 = 12+1:24;
    tidx1 = 0+1:12;

    ncycle = (tidx2(1)-tidx1(1))/12;

    gr_tt(ne) = mean(div_tt_zavg(tidx2)./div_tt_zavg(tidx1));
    gr_uu(ne) = mean(div_uu_zavg(tidx2)./div_uu_zavg(tidx1));
    gr_ww(ne) = mean(div_ww_zavg(tidx2)./div_ww_zavg(tidx1));
   
end


gr_plot = gr_tt;
gr_plot = gr_plot.^(1/ncycle);

plot_mitgcm_parm = [1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10]

%%
figure(2)
clf;
set(gcf,'Color','w')
s1 = semilogx(plot_mitgcm_parm,gr_plot,'LineWidth',2);
set(gca,'Fontsize',fontsize+2)
grid on;grid minor;
% title('Mean Floquet exponents \mu (\lambda = 400 m)','Fontsize',fontsize+5)
xlabel('White noise level (^oC)','Fontsize',fontsize+5)
ylabel('Instability Growth Rate')
% legend([p1,s1],'Floquet \mu_k^b', '\mu_k^b of MITgcm','Fontsize',fontsize+5,'Position',[0.2000 0.6655 0.1676 0.1747])
% legend([p1,s1],'Floquet \mu_k^\psi', '\mu_k^w of MITgcm','Fontsize',fontsize+5,'Position',[0.2000 0.6655 0.1676 0.1747])



