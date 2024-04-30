%%% 
%%% calc_GrowthRate_MITgcm_RMSE.m
%%% 
%%% Calculate the instability growth rate of the MITgcm simulations


clear;
for  ne = 1:11
load_all

% Ntide = 20;
% tidx = 1:Ntide*12;
% No = nDumps-1;
No =  240;
tidx = 1:No;
Nt = length(tidx);
Hshear = 250;
dz = delR(end);
Nshear = round(Hshear/dz);
zidx = Nr-Nshear:Nr-1;
% zidx = 1:Nr;
% Nshear = length(zidx)
botZ =-1500;
hab_shear = zz(zidx)-botZ;

div_tt = zeros(Nt,Nshear);
div_uu = zeros(Nt,Nshear);
div_vv = zeros(Nt,Nshear);
div_ww = zeros(Nt,Nshear);

for o = tidx
    nIter = dumpIters(o);
    time_h(o) = nIter.*deltaT./3600;
    time_tidal(o) = time_h(o)/12;
    
    % tt = squeeze(rdmds([exppath,'/results/THETA_inst'],nIter));
    % uu = squeeze(rdmds([exppath,'/results/UVEL_inst'],nIter));
    % vv = squeeze(rdmds([exppath,'/results/VVEL_inst'],nIter));
    % ww = squeeze(rdmds([exppath,'/results/WVEL_inst'],nIter));
    tt = squeeze(rdmds([exppath,'/results/THETA'],nIter));
    uu = squeeze(rdmds([exppath,'/results/UVEL'],nIter));
    % vv = squeeze(rdmds([exppath,'/results/VVEL'],nIter));
    % ww = squeeze(rdmds([exppath,'/results/WVEL'],nIter));
    tt_shear = tt(:,zidx);
    uu_shear = uu(:,zidx);
    % vv_shear = vv(:,zidx);
    % ww_shear = ww(:,zidx);
    
    mean_tt_shear = mean(tt(:,zidx));
    mean_uu_shear = mean(uu(:,zidx));
    % mean_vv_shear = mean(vv(:,zidx));
    % mean_ww_shear = mean(ww(:,zidx));

    mean_tt_shear_2d = repmat(mean_tt_shear,[Nx 1]);
    mean_uu_shear_2d = repmat(mean_uu_shear,[Nx 1]);
    % mean_vv_shear_2d = repmat(mean_vv_shear,[Nx 1]);
    % mean_ww_shear_2d = repmat(mean_ww_shear,[Nx 1]);

    % for kk = 1:Nshear
        % div_tt(o,kk) = rms(tt_shear(:,kk) - mean_tt_shear(kk));
        % div_uu(o,kk) = rms(uu_shear(:,kk) - mean_uu_shear(kk));
        % div_vv(o,kk) = rms(vv_shear(:,kk) - mean_vv_shear(kk));
        % div_ww(o,kk) = rms(ww_shear(:,kk) - mean_ww_shear(kk));
    % end
        div_tt(o,:) = rmse(tt_shear,mean_tt_shear_2d,1);
        div_uu(o,:) = rmse(uu_shear,mean_uu_shear_2d,1);
        % div_vv(o,:) = rmse(vv_shear,mean_vv_shear_2d,1);
        % div_ww(o,:) = rmse(ww_shear,mean_ww_shear_2d,1);

end

%%

CLIM = [0 1]*1e-6;

% % figure(1)
% % % clf;set(gcf,'Color','w','Position', [75 224 1362 647])
% % % subplot(2,2,1)
% % % pcolor(time_h,hab_shear,div_tt');shading flat;colorbar
% % % ylabel('HAB (m)');xlabel('Time (hours)');
% % % title('RMSE of potential temperature (degC)')
% % % set(gca,'Fontsize',fontsize)
% % % clim(CLIM)
% % subplot(2,2,2)
% % pcolor(time_h,hab_shear,div_uu');shading flat;colorbar
% % ylabel('HAB (m)');xlabel('Time (hours)');
% % set(gca,'Fontsize',fontsize)
% % title('RMSE of cross-isobath velocity u RMSE (m/s)')
% % clim(CLIM)
% % % subplot(2,2,3)
% % % pcolor(time_h,hab_shear,div_vv');shading flat;colorbar
% % % ylabel('HAB (m)');xlabel('Time (hours)');
% % % set(gca,'Fontsize',fontsize)
% % % title('RMSE of along-isobath velocity v (m/s)')
% % % clim(CLIM)
% % subplot(2,2,4)
% % pcolor(time_h,hab_shear,div_ww');shading flat;colorbar
% % ylabel('HAB (m)');xlabel('Time (hours)');
% % set(gca,'Fontsize',fontsize)
% % title('RMSE of vertical velocity w (m/s)')
% % colormap(WhiteBlueGreenYellowRed(0))
% % clim(CLIM)
% % 
% % % % % % print('-dpng','-r150',[expname '_rmse.png']);

%%
div_tt_zavg = mean(div_tt,2);
div_uu_zavg = mean(div_uu,2);
% div_vv_zavg = mean(div_vv,2);
% div_ww_zavg = mean(div_ww,2);

div_tt_norm = div_tt_zavg/div_tt_zavg(1);
div_uu_norm = div_uu_zavg/div_uu_zavg(1);
% div_vv_norm = div_vv_zavg/div_vv_zavg(1);
% div_ww_norm = div_ww_zavg/div_ww_zavg(1);




%%% Calculate the growth rate
    fit_span = 12*1+1:10*12;
    xxplot = time_h;
    yyplot = log(div_uu_zavg/2)/2;
    [pp,S] = polyfit(xxplot(fit_span),yyplot(fit_span),1); 
    grow = pp(1)
    [y_fit,delta_fit] = polyval(pp,xxplot,S);

figure()
clf;set(gcf,'Color','w','Position',[211 289 852 394])
plot(time_h/12,log(div_tt_zavg/2)/2,'LineWidth',2);
hold on;
plot(time_h/12,log(div_uu_zavg/2)/2,'LineWidth',2);
% plot(xxplot(fit_span)/12, y_fit(fit_span),'--');
set(gca,'Fontsize',fontsize)
xlabel('Time (tidal cycles)')
title('Normalized temperature RMSE')
title('RMSE of T and u averaged over the bottom shear layer')
ylabel('(degC)')
grid on;grid minor;
hold on;
% ylim([1e-9 1e-1])
% ylim([min(min([div_tt_zavg div_uu_zavg])) max(max([div_tt_zavg div_uu_zavg]))])

% plot(div_uu_norm);
% plot(div_vv_norm);
% plot(div_ww_norm);
% % % % print('-dpng','-r150',[expname '_rmse2.png']);


% filename = [expdir expname '/RMSE_mean.mat'];

% save(filename,'time_h','hab_shear','div_tt_zavg','div_uu_zavg','div_vv_zavg','div_ww_zavg',...
%     'div_tt_norm','div_vv_norm','div_ww_norm','div_uu_norm',...
%     'div_tt','div_uu','div_vv','div_ww')

% save(filename,'time_h','hab_shear','div_tt_zavg','div_uu_zavg','div_ww_zavg',...
%     'div_tt_norm','div_ww_norm','div_uu_norm',...
%     'div_tt','div_uu','div_ww')


end

