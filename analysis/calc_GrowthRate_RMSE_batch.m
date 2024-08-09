%%% 
%%% calc_GrowthRate_RMSE.m
%%% 
%%% Calculate the instability growth rate of the MITgcm simulations


clear;
close all



for  ne = 1:22

    NTtide = NTtide_all(ne);
    Nstart = Nstart_all(ne);
    Nend = Nend_all(ne);

    load_all;
    No = NTtide*12;
    tidx = 1:No;
    
    fit_span = Nstart*12+1:Nend*12;
    
    Nt = length(tidx);
    dz = delR(end);
    zidx = 351:500;
    Nshear = length(zidx);
    
    div_tt = zeros(Nt,Nshear);
    div_ww = zeros(Nt,Nshear);
    
    parfor o = tidx
        nIter = dumpIters(o);
        time_h(o) = nIter.*deltaT./3600;
        time_tidal(o) = time_h(o)/12;
    
        tt = squeeze(rdmds([exppath,'/results/THETA'],nIter));
        tt_shear = tt(:,zidx);
        mean_tt_shear = mean(tt(:,zidx));
        mean_tt_shear_2d = repmat(mean_tt_shear,[Nx 1]);
        div_tt(o,:) = rmse(tt_shear,mean_tt_shear_2d,1);
    
        ww = squeeze(rdmds([exppath,'/results/WVEL'],nIter));
        ww_shear = ww(:,zidx);
        mean_ww_shear = mean(ww(:,zidx));
        mean_ww_shear_2d = repmat(mean_ww_shear,[Nx 1]);
        div_ww(o,:) = rmse(ww_shear,mean_ww_shear_2d,1);
    
    end
    
    %%
    div_tt2_zint = sum(dz*div_tt.^2,2);
    div_ww2_zint = sum(dz*div_ww.^2,2);
    
    Pr = 1;
    ke = div_ww2_zint/2;
    pe = Pr*div_tt2_zint/2;
    
    xxplot = time_h;
    yyplot = log(pe)/2;
    [pp,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
    grow(ne) = pp(1)
    [y_fit,delta_fit] = polyval(pp,xxplot,S);
    
    
    filename = [expdir expname '/RMSE_tt.mat'];
    save(filename,'time_h','xxplot','yyplot','fit_span','pp','y_fit',...
        'div_ww2_zint','div_tt2_zint','grow')
    
    
    
    h_figure = figure(1);
    clf;
    set(gcf,'Color','w','Position',[211 289 852 394])
    plot(time_h/12,log(pe)/2,'LineWidth',2);
    hold on;
    plot(time_h/12,log(ke)/2,'LineWidth',2);
    plot(xxplot(fit_span)/12, y_fit(fit_span),'k-.','LineWidth',2);
    set(gca,'Fontsize',fontsize+5)
    xlabel('Time (tidal cycles)','Interpreter','latex')
    title('Normalized turbulent energy in the shear layer','Interpreter','latex')
    ylabel('$\log(\mathrm{turbulent\ energy})$','Interpreter','latex')
    legend('Turbulent potential energy','Turbulent kinetic energy',...
        'Linear fit','Position',[0.5802 0.2367 0.3174 0.1832]...
        ,'Interpreter','latex','Fontsize',fontsize-4)
    grid on;grid minor;
    hold on;
    print('-dpng','-r150',[expdir expname '_rmse_tt.png']);


end
